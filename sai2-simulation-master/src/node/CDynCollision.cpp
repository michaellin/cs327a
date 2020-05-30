//===========================================================================

/*
    This file is part of the dynamics library.
    Copyright (C) 2014, Artificial Intelligence Laboratory,
    Stanford University. All rights reserved.

    \version   $MAJOR.$MINOR.$RELEASE $Rev: 1242 $
*/
//===========================================================================

//---------------------------------------------------------------------------
#include "utility/CDynLogger.h"
#include "object/CDynObject.h"
#include "node/CDynWorld.h"
#include "node/CDynContact.h"
#include "node/CDynCollisionCheckRecord.h"
#include "node/CDynCollision.h"
#include "distance/CDynBSphere.h"
//---------------------------------------------------------------------------
cDynPair* cDynPair::free_=NULL;
static const cDynTime tlimit=1e-8;
//---------------------------------------------------------------------------
int g_dynPairCheck=0;
//---------------------------------------------------------------------------
//extern bool g_dynNoBackup;
//---------------------------------------------------------------------------

cDynCollision::cDynCollision()
{
    {for (int i=0;i<3;i++) 
    {
        head_[i]=tail_[i]=NULL;
    }}
    {for (int i=0;i<DEPAIR_HASHSIZE;i++) 
    {
        hash_[i]=NULL;
    }}
    {for (int i=0;i<CDYN_PAIR_NUMSTATES;i++) 
    {
        num_[i]=0;
        list_[i]=NULL;
    }}
    reset_=true;
}

//---------------------------------------------------------------------------

cDynCollision::~cDynCollision()
{
    // {for (int i=0;i<3;i++) 
    // {
    //     assert(head_[i] == NULL);
    //     assert(tail_[i] == NULL);
    // }}
    // {for (int i=0;i<DEPAIR_HASHSIZE;i++) 
    // {
    //     assert(hash_[i] == NULL);
    // }}
    // {for (int i=0;i<CDYN_PAIR_NUMSTATES;i++) 
    // {
    //     assert(num_[i] == 0);
    //     assert(list_[i] == NULL);
    // }}
}

//---------------------------------------------------------------------------

void cDynCollision::add(cDynObject* obj)
{
    assert(obj->geometry.box() == NULL);
    assert(obj->baseNode()->world() != NULL);
    const cDynBSphere* bs=obj->geometry.bs();
    if (bs == NULL) return;

    cDynTime time=obj->baseNode()->world()->time();

    cDynBoundRecord* b = new cDynBoundRecord;
    obj->geometry.box(b);

    double range[2][3];
    if (obj->isFixed()) 
    {
        cDynTransform T;
        T.set(obj->globalFrame(time));
        cDynVector3 col;
        cDynPrim* prim;
        int w;
        for (int axis=0;axis<3;axis++) 
        {
            prim=NULL;
            col.set(-T.rotation(axis,0),-T.rotation(axis,1),-T.rotation(axis,2));
            range[0][axis]= -bs->max(&col,prim, &w,0,NULL) 
                            + T.translation(axis) - (double)1e-6;
            prim=NULL;
            col.negate(col);
            range[1][axis]= bs->max(&col,prim, &w,0,NULL)
                            + T.translation(axis) + (double)1e-6;
        }
    } 
    else 
    {
        cDynVector3 p;
        p.multiply(obj->globalFrame(time), CDYN_BSPHERE_P(bs));
        double r=CDYN_BSPHERE_R(bs);
        for (int i=0;i<3;i++) 
        {
            range[0][i]=p[i]-r;
            range[1][i]=p[i]+r;
        }
    }

    {for (int axis=0;axis<3;axis++) 
    {
        // set values
        cDynBound* min=&(b->min_[axis]);
        cDynBound* max=&(b->max_[axis]);
        min->type_= -1;
        max->type_=  1;

        min->obj_= obj;
        max->obj_= obj;

        min->set(range[0][axis],time);
        max->set(range[1][axis],time);

        if (head_[axis] == NULL) 
        {
            // first bound
            min->left_=NULL;
            min->right_=max;
            max->left_=min;
            max->right_=NULL;
            head_[axis]=min;
            tail_[axis]=max;
        } 
        else 
        {
            // find position in list for min bound
            cDynBound* pos=head_[axis];
            cDynBound* prev=NULL;
            while (pos != NULL && pos->value_ < min->value_) 
            {
                pairUpdate(obj,pos->obj_, -pos->type_);
                prev=pos;
                pos=pos->right_;
            }

            // insert min bound into list. Update head & tail if needed.
            min->left_=prev;
            min->right_=pos;
            if (prev == NULL) head_[axis]=min;
            else		  prev->right_=min;
            
            if (pos == NULL)  tail_[axis]=min;
            else		  pos->left_=min;

            // insert max bound
            insert(axis, max, min, pos, 1);
        }
    }}
    //check();
}

//---------------------------------------------------------------------------

void cDynCollision::update(cDynObject* node, const cDynTime& min, const cDynTime& max)
{
    cDynBoundRecord* b=node->geometry.box();
    if (b == NULL) return;
    //check();
    const cDynBSphere* bs=node->geometry.bs();
    assert(bs != NULL);
    cDynVector3 p;
    p.multiply(node->globalFrame(max), CDYN_BSPHERE_P(bs));
    double r=CDYN_BSPHERE_R(bs);

    for (int axis=0;axis<3;axis++) 
    {
        for (int i=0;i<2;i++) 
        {
            cDynBound* bound;
            if (i==0) 
            {
                bound=&(b->min_[axis]);
                bound->update(p[axis]-r,min,max);
            } 
            else 
            {
                bound=&(b->max_[axis]);
                bound->update(p[axis]+r,min,max);
            }

            if (bound->left_ != NULL && bound->left_->value_ > bound->value_) 
            { 
                // bound shifts down
                // remove from list
                bound->left_->right_=bound->right_;
                if (bound->right_ != NULL)
                    bound->right_->left_=bound->left_;
                else
                    tail_[axis]=bound->left_;
                if (bound->type_ != bound->left_->type_)
                    pairUpdate(bound->obj_,bound->left_->obj_,bound->left_->type_);

                // reinsert
                insert(axis, bound, bound->left_->left_, bound->left_, -1);
            } 
            else if (bound->right_ != NULL && bound->right_->value_ < bound->value_) 
            { 
                // bound shifts up
                // remove from list
                bound->right_->left_=bound->left_;
                if (bound->left_ != NULL)
                    bound->left_->right_=bound->right_;
                else
                    head_[axis]=bound->right_;
                if (bound->type_ != bound->right_->type_)
                    pairUpdate(bound->obj_,bound->right_->obj_,bound->type_);
                // reinsert
                insert(axis, bound, bound->right_, bound->right_->right_, 1);
            }
        }
    }
}

//---------------------------------------------------------------------------

void cDynCollision::insert(const int axis, cDynBound* bound, cDynBound* prev, cDynBound* next, const int dir)
{
    if (dir == 1) 
    {
        while (next != NULL && next->value_ < bound->value_) 
        {
            if (bound->type_ != next->type_) pairUpdate(bound->obj_, next->obj_, bound->type_);
            prev=next;
            next=next->right_;
        }
    } 
    else 
    { 
        // dir == -1
        while (prev != NULL && prev->value_ > bound->value_) 
        {
            if (bound->type_ != prev->type_) pairUpdate(bound->obj_, prev->obj_, prev->type_);
            next=prev;
            prev=prev->left_;
        }
    }

    bound->left_=prev;
    bound->right_=next;

    if (prev == NULL) head_[axis]=bound;
    else             prev->right_=bound;
    if (next == NULL) tail_[axis]=bound;
    else		  next->left_=bound;
}

//---------------------------------------------------------------------------

cDynPair* cDynCollision::pairUpdate(cDynObject* a, cDynObject* b, const int inc, const int s)
{
    if (a == b) return(NULL);

    unsigned long al=a->uid()+1;
    unsigned long bl=b->uid()+1;

    unsigned long h=(al * bl)%DEPAIR_HASHSIZE;

    cDynPair* p=hash_[h];
    cDynPair* o=NULL;

    while (p != NULL) {
        if ((p->a_ == a && p->b_ == b)||(p->b_ == a && p->a_ == b)) 
        {
            if (s < 0) return(p);
            if (inc == 0) 
            {
                if (s == 0) state(p,(p->count_ == 3)?CDYN_PAIR_ACTIVE:CDYN_PAIR_INACTIVE);
                else state(p, s);
            }

            p->count_ += inc;
            if (p->count_ == 0 && (p->state_ == CDYN_PAIR_ACTIVE || p->state_ == CDYN_PAIR_INACTIVE)) 
            {
                // if zero remove from hash
                if (o == NULL) hash_[h]=p->hash_;
                else o->hash_=p->hash_;
                state(p, -1);
                delete p;
                p=NULL;
            } 
            
            else if (p->state_ == CDYN_PAIR_ACTIVE && p->count_ != 3) 
            {
                // pair no longer in active list
                state(p, CDYN_PAIR_INACTIVE);
            } 
            
            else if (p->state_ == CDYN_PAIR_INACTIVE && p->count_ == 3) 
            {
                state(p, CDYN_PAIR_ACTIVE);
            }
            return(p);
        }
        o=p;
        p=p->hash_;
    }

    if (s < 0) return(NULL);
    
    // create a new pair to 
    if (a->isFixed() && b->isFixed()) return(NULL);
    p=new cDynPair;
    
    if (a->isFixed()) 
    { 
        // if one is fixed it will always be stored in b_
        p->a_=b;
        p->b_=a;
    } 
    else 
    {
        p->a_=a;
        p->b_=b;
    }
    p->count_=inc;

    // add to state list
    if (s == 0) p->state_=(p->count_ == 3)?CDYN_PAIR_ACTIVE:CDYN_PAIR_INACTIVE;
    else p->state_=s;

    // add it to the appropriate state list
    p->prev_=NULL;
    p->next_=list_[p->state_];
    if (p->next_ != NULL) p->next_->prev_=p;
    list_[p->state_]=p;
    num_[p->state_]++;

    // add to hash table
    p->hash_=hash_[h];
    hash_[h]=p;

    return(p);
}

//---------------------------------------------------------------------------

void cDynCollision::ignore(cDynObject* a, cDynObject* b, const bool cond)
{
    if (cond)
        pairUpdate(a,b,0,CDYN_PAIR_IGNORE);
    else
        pairUpdate(a,b,0,0);
}

//---------------------------------------------------------------------------

void cDynCollision::invalid(cDynObject* a, cDynObject* b, const bool cond)
{
    if (cond)
        pairUpdate(a,b,0,CDYN_PAIR_INVALID);
    else
        pairUpdate(a,b,0,0);
}

//---------------------------------------------------------------------------

void cDynCollision::activate(cDynObject* a, cDynObject* b, const bool cond)
{
    if (cond)
        pairUpdate(a,b,0,0);
    else
        pairUpdate(a,b,0,CDYN_PAIR_SLEEP);
}

//---------------------------------------------------------------------------

int cDynCollision::pairState(cDynObject* a, cDynObject* b)
{
    cDynPair* p = pairUpdate(a,b,0,CDYN_PAIR_QUERY);
    return((p == NULL)?CDYN_PAIR_INACTIVE:p->state_);
}

//---------------------------------------------------------------------------

void cDynCollision::state(cDynPair* p, const int state)
{
    if (p->state_ == state) return;

    // disconnect from current list
    if (p->prev_ == NULL) 
    {
        list_[p->state_]=p->next_;
        if (list_[p->state_] != NULL) list_[p->state_]->prev_=NULL;
    } 
    else 
    {
        p->prev_->next_=p->next_;
        if (p->next_ != NULL) p->next_->prev_=p->prev_;
    }

    num_[p->state_]--;
    // change state
    p->state_=state;
    if (state < 0) return; // remove from list

    if (state == CDYN_PAIR_ACTIVE && p->a_->baseNode()) 
    {
        if (p->safe_ < p->a_->baseNode()->world()->time())
            p->safe_ = p->a_->baseNode()->world()->time(); // XXX check if this is the best time
    }

    // add to new list
    p->prev_=NULL;
    p->next_=list_[p->state_];
    if (list_[p->state_] != NULL) list_[p->state_]->prev_=p;
    list_[p->state_]=p;
    num_[p->state_]++;
    assert(p != p->next_);
    assert(p != p->prev_);
}

//---------------------------------------------------------------------------

void cDynCollision::remove(cDynObject* node)
{
    cDynBoundRecord* r=node->geometry.box();
    if (r != NULL) 
    {
        for (int axis=0;axis<3;axis++) 
        {
            cDynBound* min=&(r->min_[axis]);
            cDynBound* max=&(r->max_[axis]);
            assert(max->left_ != NULL); // min has to be to the left of max
            
            // disconnect max from list
            max->left_->right_=max->right_;
            if (max->right_ == NULL) tail_[axis]=max->left_;
            else  max->right_->left_=max->left_;
            
            // shift till you hit min
            {for (cDynBound* b=max->left_;b != min; b=b->left_) 
            {
                if (max->type_ != b->type_)
                    pairUpdate(node, b->obj_, b->type_);
            }}

            // disconnect min from list
            if (min->left_  == NULL) head_[axis]=min->right_;
            else min->left_->right_=min->right_;
            if (min->right_ == NULL) tail_[axis]=min->left_;
            else min->right_->left_=min->left_;
            for (cDynBound* b=min->left_;b != NULL; b=b->left_)
                pairUpdate(node, b->obj_, b->type_);
        }
        delete r;
        node->geometry.box(NULL);
    }
}

//---------------------------------------------------------------------------

void cDynCollision::pop(cDynPair* p)
{
    // pop pair to top of state list

    // disconnect from list list
    if (p->prev_ == NULL) 
    {
        return; // already at top of list
    } 
    else
    {
        p->prev_->next_=p->next_;
        if (p->next_ != NULL) p->next_->prev_=p->prev_;
    }

    // add at top
    p->prev_=NULL;
    p->next_=list_[p->state_];
    if (list_[p->state_] != NULL) list_[p->state_]->prev_=p;
    list_[p->state_]=p;
    assert(p != p->next_);
    assert(p != p->prev_);
}

//---------------------------------------------------------------------------

int cDynCollision::checkObjects(const cDynTime& time, cDynTestOption option, cDynBaseNode* a)
{
    cDynPair* p=list_[CDYN_PAIR_ACTIVE];
    while (p != NULL) 
    {
        cDynPair* l=p;
        p=p->next_; // update early in case options change list
        if (a == NULL || l->a_->baseNode() == a || l->b_->baseNode() == a) 
        {
            cDynPrimPairCallnum();
            cDynCollisionCheckRecord rec(CDYN_ALL_CONTACT);
            cDynFrame F;
            F.inversedMultiply(l->a_->globalFrame(time), l->b_->globalFrame(time));
            cDynTransform Tab;
            Tab.set(F);
            int type=rec.checkBS2(l->a_,l->b_,Tab);
            if (type == 0) 
            {
                if (option == CDYN_TESTOPTION_PRINT) 
                {
                    cDynPrintf("    contact %30s <--> %-30s\n",
                        cDynObjectName(l->a_),
                        cDynObjectName(l->b_)
                    );
                } 
                else if (option == CDYN_TESTOPTION_IGNORE) 
                { 
                    // added because objects may be in contact but penetrating due to error term
                    cDynPrimPairCallnum();
                    cDynPrimPair* pair=rec.distance(l->a_,l->b_,Tab);
                    if (pair->distance() < 1e-10) 
                    {
/*
                        cDynPrintf("ignoring contact at %30s <--> %-30s\n",

                        cDynObjectName(l->a_),
                        cDynObjectName(l->b_)
                        );
*/
                        //ignore(l->a_, l->b_);
                    }
                }
            } 
            else if (type == -1) 
            {
                if (option == CDYN_TESTOPTION_VALID) 
                {
                    return(-1);
                } 
                else if (option == CDYN_TESTOPTION_IGNORE) 
                {
    /*
                    cDynPrintf("ignoring %30s <--> %-30s\n",

                        cDynObjectName(l->a_),
                        cDynObjectName(l->b_)
                    );
    */
                    ignore(l->a_, l->b_);
                } 
                else if (option == CDYN_TESTOPTION_PRINT) 
                {
                    cDynPrintf("PENETRATION %30s <--> %-30s\n",
                        cDynObjectName(l->a_),
                        cDynObjectName(l->b_)
                    );
                }
            }
        }
    }
    return(1);
}

//---------------------------------------------------------------------------

void cDynCollision::checkInvalids()
{
    cDynPair* p=list_[CDYN_PAIR_INVALID];
    bool jump=false;
    while (p != NULL) 
    {
        cDynTime t=p->a_->baseNode()->world()->time();

        cDynPrimPairCallnum();
        cDynCollisionCheckRecord rec(CDYN_ALL_CONTACT);
        cDynFrame F;
        F.inversedMultiply(p->a_->globalFrame(t), p->b_->globalFrame(t));
        cDynTransform Tab;
        Tab.set(F);
        int type=rec.checkBS2(p->a_,p->b_,Tab);

        if (type == 1) 
        { 
            // type == 1 means no penetration return it to active/inactive list
            cDynPair* l=p; p=p->next_;
            jump=true;
            cDynPrintf("reactivating collisions between %s and %s\n",cDynObjectName(l->a_),cDynObjectName(l->b_));
            invalid(l->a_,l->b_,false);
        }

        if (jump) jump=false;
        else p=p->next_; 
    }
}

//---------------------------------------------------------------------------

void cDynCollision::checkSleep(const cDynBaseNode* b)
{
    for (cDynPair*p=list_[CDYN_PAIR_ACTIVE];p != NULL;p=p->next_) 
    {
        if (p->a_->baseNode() == b || p->b_->baseNode() == b) 
        {
            p->safe_ = b->world()->time();
        }
    }
}

//---------------------------------------------------------------------------

bool cDynCollision::checkContacts(cDynPair* start, cDynPair* end, cDynBaseNode* a, cDynBaseNode* b)
{
    static int flag=0;
    static int iter=0;

    flag=(cDynDebugLevel > 0)?1:0;
    iter++;

    bool all=(a == NULL && b == NULL);
    if (start == NULL) start=list_[CDYN_PAIR_ACTIVE];
    //check();
    if (flag) 
    {
        cDynPrintf("COLLISION %d\n",iter);
    }

    cDynPair* n=NULL;
    for (cDynPair* p=start; p != end; p=n) 
    {
        n=p->next_;
        if (all 
         || p->a_->baseNode() == a || p->a_->baseNode() == b 
         || p->b_->baseNode() == a || p->b_->baseNode() == b ) 
        {
            cDynTime tmax,tA,tB;
            if (p->a_->baseNode()->status() == CDYN_SLEEP) 
            {
                if (p->b_->isFixed() || p->b_->baseNode()->status() == CDYN_SLEEP) continue;
                p->flip(); // always make the sleeping node B
            } 

            if (p->b_->isFixed() || p->b_->baseNode()->status() == CDYN_SLEEP) 
            {
                tmax=tA=tB=p->a_->baseNode()->contact()->time();
            } 
            else 
            {
                tA=p->a_->baseNode()->contact()->time();
                tB=p->b_->baseNode()->contact()->time();
                tmax=(tA < tB)?tA:tB;
            }

            //if (!reset_ && tmax == p->safe_) continue;
            cDynTime tmin=(p->safe_ < tmax)?p->safe_:p->a_->baseNode()->world()->time();

            cDynPrimPairCallnum();
            cDynCollisionCheckRecord rec(CDYN_ALL_CONTACT + CDYN_SAVE_CONTACT + CDYN_SAVE_PENETRATION);
            g_dynPairCheck++;
            cDynFrame F;
            F.inversedMultiply(p->a_->globalFrame(tmax), p->b_->globalFrame(tmax));
            cDynTransform Tab;
            Tab.set(F);
            int type=rec.checkBS2(p->a_,p->b_,Tab);
            if (flag) 
            {
                //p->a_->baseNode()->displaySmall(tmax);
                cDynPrintf("%2d time=%-13.9f %s<--->%s\n", type, tmax,
                    cDynObjectName(p->a_),cDynObjectName(p->b_));
            }

            switch (type) 
            {
                case  1: p->safe_=tmax; break; // no contact, only need to update safe time
                case  0: p->safe_=tmax;
                    if (tA != tmax) 
                    {
                        p->a_->baseNode()->contact()->flush();
                        p->a_->baseNode()->contact()->updatePersistent(tmax);
                        p->a_->baseNode()->checkJointLimits(tmax);
                    }
                    if (tB != tmax) 
                    {
                        p->b_->baseNode()->contact()->flush();
                        p->b_->baseNode()->contact()->updatePersistent(tmax);
                        p->b_->baseNode()->checkJointLimits(tmax);
                    }

                    // contact exists add to list
                    rec.addContacts(tmax); 

                    // next, if the contact time changed need to recheck all pairs where time has changed
                    if (tA != tmax || tB != tmax) 
                    { 
                        // need to go over list again
                        pop(p);
                        n=p;
                    }
                    break;
                case -1:  // penetration find new contact time
#if 0
                    if (g_dynNoBackup) 
                    {
                            cDynPrintf("backup invalid %s <--> %s\n",cDynObjectName(p->a_), cDynObjectName(p->b_));
                            invalid(p->a_, p->b_);
                    } 
                    else 
                    {
#endif
                        cDynPrimPairError error=CDYN_NO_FAILURE;
                        for (cDynPrimPair* pair=rec.plist();pair != NULL && error == CDYN_NO_FAILURE; pair=pair->next()) 
                        {
                            tmax=pair->checkPrimTime(tmin,tmax,Tab,error);
                            switch (error) 
                            {
                                case CDYN_NO_FAILURE: 
                                    break;
                                
                                case CDYN_UNKNOWN_FAILURE: 
                                    cDynPrintf("UNKNOWN ITERATION FAILURE between %s and %s (%5.13f, %5.13f)\n", 
                                    cDynObjectName(p->a_), cDynObjectName(p->b_),tmin,tmax);
                                    break;

                                case CDYN_REGRESSION_FAILURE: 
                                    cDynPrintf("REGRESSION FAILURE between %s and %s\n", 
                                    cDynObjectName(p->a_), cDynObjectName(p->b_));
                                    break;

                                case CDYN_TOLERANCE_FAILURE: 
                                    cDynPrintf("TOLERANCE FAILURE between %s and %s\n", 
                                    cDynObjectName(p->a_), cDynObjectName(p->b_));
                                    break;
                            }
                        }
                        if (error != CDYN_NO_FAILURE) 
                        { 
                            // mark as invalid and check again
                            invalid(p->a_, p->b_);
                            return(false);
                        } 
                        else 
                        {
                            // move to top of list
                            pop(p);
                            p->a_->baseNode()->contact()->flush();
                            p->a_->baseNode()->contact()->updatePersistent(tmax);
                            p->a_->baseNode()->checkJointLimits(tmax);
                            if (p->b_->baseNode() != p->a_->baseNode()) 
                            {
                                if (!(p->b_->isFixed() || p->b_->baseNode()->status() == CDYN_SLEEP)) 
                                {
                                    p->b_->baseNode()->contact()->flush();
                                    p->b_->baseNode()->contact()->updatePersistent(tmax);
                                    tmax=p->b_->baseNode()->checkJointLimits(tmax);
                                }
                            }
                            n=p;
                        }
#if 0
                    }
#endif
                    break;
            } // switch
        } // if
    } // for
    reset_=false;
    return(true);
}

//---------------------------------------------------------------------------
#ifdef CDYN_DEBUG
//---------------------------------------------------------------------------

void cDynCollision::display()
{
    for (cDynPair* l=list_[CDYN_PAIR_ACTIVE];l != NULL; l=l->next_) 
    {
        cDynPrintf("%30s <--> %-30s\n",
            cDynObjectName(l->a_),
            cDynObjectName(l->b_)
        );
    }
}

//---------------------------------------------------------------------------

void cDynCollision::check()
{
    int fault=0;
    static int iter=0;
    iter++;

    {for (int i=0;i<3;i++) 
    {
        cDynPair* p=list_[i];
        while (p != NULL) 
        {
            if (p->state_ != i) fault++;
            if (p->prev_==NULL) 
            {
                if (p != list_[i]) fault++;
            } 
            else 
            {
                if (p->prev_ == p) fault++;
                if (p->prev_->next_ != p) fault++;
            }	
            if (p->next_!=NULL) 
            {
                if (p->next_ == p) 
                {
                    cDynPrintf("iter=%d\n",iter);
                    exit(1);
                }
                if (p->next_->prev_ != p) fault++;
            }	
                
            p=p->next_; // update early in case options change list
        }
    }}
    {for (int i=0;i<DEPAIR_HASHSIZE;i++) 
    {
        cDynPair* c=hash_[i];
        while (c != NULL) 
        {
            cDynPair* p=list_[c->state_];
            while (p != NULL) 
            {
                if (p == c) break;
                p=p->next_;
            }
            if (p != c) fault++;
            c=c->next_;
        }
    }}
    {for (int axis=0;axis<3;axis++) 
    {
        cDynBound* b=head_[axis];
        while (b != NULL) 
        {
            if (b==b->left_) fault++;
            if (b==b->right_) fault++;
            if (b->left_ == NULL) 
            {
                if (b!=head_[axis]) fault++;
            } 
            else 
            {
                if (b->left_->value_ > b->value_) fault++;
                if (b->left_->right_ != b) fault++;
            }

            if (b->right_ == NULL) 
            {
                assert(b==tail_[axis]);
            } 
            else 
            {
                if (b->right_->value_ < b->value_) fault++;
                if (b->right_->left_ != b) fault++;
            }
            b=b->right_;
        }
    }}
    if (fault) 
    {
        cDynPrintf("iter=%d\n",iter);
        exit(1);
    }
}

//---------------------------------------------------------------------------
#endif // CDYN_DEBUG
//---------------------------------------------------------------------------

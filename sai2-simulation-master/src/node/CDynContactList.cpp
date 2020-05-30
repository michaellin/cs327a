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
//---------------------------------------------------------------------------
// XXX used for rand() during testing
#undef DEBUG
#include <math.h>
#include <stdlib.h>
#include <string.h>
#ifdef BORLAND
#include <mem.h>
#endif
#include "matrix/CDynMathDefn.h"
#include "node/CDynContact.h"
#include "node/CDynContactEvent.h"
#include "node/CDynConstraint.h"
#include "dynamics/CDynDynamics.h"
#include "node/CDynWorld.h"
#include "object/CDynObject.h"
#include "node/CDynBaseNode.h"
//---------------------------------------------------------------------------
#undef RTIMER1
#ifdef RTIMER1
#include "utility/CDynTimeInterval.h"
extern cDynTimeInterval rtimer;
extern cDynTime rtotal;
extern int rcount;
#endif
//---------------------------------------------------------------------------
#undef CHECKMATH
//---------------------------------------------------------------------------
#ifdef C_DYN_DEBUG_OPS
int const ESIZE=2;
double KE[ESIZE];
double PE[ESIZE];
double DeltaEnergy;
#endif
//---------------------------------------------------------------------------
#define MAXSTACKALLOC	(1e6)
//---------------------------------------------------------------------------
//#include "matrix/CDynMathDefn.h" // CDYN_ERROR_BOUND
#define CDYN_ERROR_BOUND 1e-6f
#define VLIMIT  1e-4
//---------------------------------------------------------------------------
static const int CDYN_UNILATERAL = 0;
static const int CDYN_BILATERAL  = 1;
static const int CDYN_REDUNDANT  = 2;
static const int CDYN_NOTACTIVE  = 3;
//---------------------------------------------------------------------------
static const int CDYN_FRICTION_NORMAL    = 4;
static const int CDYN_SATURATED_NORMAL    = 8;
static const int CDYN_FRICTION   = 16;
//---------------------------------------------------------------------------
#define BASE(x) ((x)&0x03)
#define FLAG(x) ((x)&(CDYN_FRICTION+CDYN_FRICTION_NORMAL+CDYN_SATURATED_NORMAL))
//---------------------------------------------------------------------------
//enumeration of passes in RESOLVE
static const int CDYN_IMPULSE = 0;
static const int CDYN_FORCE = 1;
static const int CDYN_ERROR = 2;
static const int CDYN_DONE = 3;
//---------------------------------------------------------------------------
static const int CDYN_NOTSET = -1;
//---------------------------------------------------------------------------
static int currentSize;
static int currentNum;
static int currentLambdaNum;
static int currentCol;
static cDynContact** currentContact;
static cDynContactPoint* currentContactPoint;
static double currentContactPointDirection=0.0f;
static double* currentLambda;
static double* currentX;
static cDynVector3 retV;	// return vector used by cDynCallbackJtSphereIn
static cDynVector6 currentS;
static int*	currentIndex;
//---------------------------------------------------------------------------
// used during debugging
static double* testX;
static double* testA;
static double* testB;
static double* testE;
static double* testErr;
static double* testV;
static double* testDV;
static double* testDelta;
static double* testDiag;
static double* testMax;
//static (int (*)[3])  testType;
static int*	testIndex;
static int pass;
//---------------------------------------------------------------------------

static inline double& lambda(const int i, const int j) 
{
    return((i<j)?
         currentLambda[i+currentSize*j]
        :currentLambda[j+currentSize*i]
    );
}

//---------------------------------------------------------------------------

// lower half of Lambda including diag.

static inline double& lcol(const int j) 
{ 
  return(currentLambda[currentCol+j]); 
}

//---------------------------------------------------------------------------

static inline double& lambda(const int i, const int j, const int* index) 
{
    return((index[i] < index[j])?	//only access upper half of matrixm --> lower half + diag
         currentLambda[index[i]+currentSize*index[j]]
        :currentLambda[index[j]+currentSize*index[i]]
    );
}
#if 0
static inline void prefetchlambda(const int i, const int j, const int* index) 
{
    prefetch((index[i] < index[j])?	//only access upper half of matrixm --> lower half + diag
         &currentLambda[index[i]+currentSize*index[j]]
        :&currentLambda[index[j]+currentSize*index[i]]
    );
}
#endif

static inline double& L(const int i, const int j, const int* index)
{
    return((index[i] < index[j])?	// only access lower half, diag not valid --> upper half
         currentLambda[index[j]+currentSize*index[i]]
        :currentLambda[index[i]+currentSize*index[j]]
    );
};

//---------------------------------------------------------------------------
#ifdef C_DYN_DEBUG_OPS
//---------------------------------------------------------------------------

void cDynDebugEnergyReset()
{
    for (int i=0;i<ESIZE;i++) 
    {
        PE[i]=0.0f;
        KE[i]=0.0f;
    }
    DeltaEnergy=0.0f;
}

//---------------------------------------------------------------------------

void cDynDebugEnergy(cDynBaseNode* b, int i)
{
    b->dynamics()->fwdDynamics();
    PE[i] += b->dynamics()->potentialEnergy();
    KE[i] += b->dynamics()->kineticEnergy();
}

//---------------------------------------------------------------------------

double cDynDebugKE(int i) { return(KE[i]); }
double cDynDebugPE(int i) { return(PE[i]); }
double cDynDebugTotal(int i) { return(PE[i]+KE[i]); }
double cDynDebugDeltaEnergy(double e) {return(((e*e)-1.0f)*DeltaEnergy);}

//---------------------------------------------------------------------------
#endif
//---------------------------------------------------------------------------

cDynContactList::cDynContactList(cDynBaseNode* base, const int size)
{
    base_=base;
    //	list_=new (cDynContact *)[size];
    list_=new cDynContact*[size];
//	if (cDynDebugLevel == -100)
//		cDynPrintf("cDynContactList\n");

    for (int i=0;i<size;i++)
        list_[i] = NULL;
    contactEvent_=NULL;
    size_=size;
    time_=CDYN_TIME_MAX;
    error_=false;
    wakeForce_=1e-8;
    wakeImpulse_=1e-8;
    n_=0;

    setup_=CDYN_SETUP_UNINITIALIZED;
    inc_ = -100;
}

//---------------------------------------------------------------------------

cDynContactList::~cDynContactList()
{
    reset();
    delete [] list_;
}

//---------------------------------------------------------------------------

void cDynContactList::grow(const int size)
{
    cDynContact** list;
    //	list=new (cDynContact *)[size];
    list=new cDynContact*[size];
//	if (cDynDebugLevel == -100)
//		cDynPrintf("cDynContactList:grow\n");

    {for (int i=0;i<size;i++)
        list[i]=NULL;}
    {for (int i=0;i<n_;i++)
        list[i]=list_[i];}
    size_=size;
    delete [] list_;
    list_=list;
}

//---------------------------------------------------------------------------

void cDynContactList::flush()
{
    //cDynPrintf("flush (%d)\n",n_);
    int n=n_;
    n_=0;
    for (int i=0;i<n;i++) 
    {
        cDynContact* c=list_[i];
        if (c->type_ == CDYN_CONTACT && c->pair_ != NULL && c->pair_->list_->n_ > 0) 
        {
            c->pair_->pair_=NULL;
            c->pair_->list_->flush();
        }

        if (c->type_ == CDYN_CONTACT && c->prim_ != NULL) 
        {
            delete c->prim_;
            c->prim_=NULL;
        }

//		list_[i]=NULL;
        delete c;
    }
    //time_=CDYN_TIME_MAX;
}

//---------------------------------------------------------------------------

void cDynContactList::reset()
{
    // flush the contact list
    for (int i=0;i<n_;i++) 
    {
        remove(i);
    }
    time_=CDYN_TIME_MAX;
    n_=0;
}

//---------------------------------------------------------------------------

void cDynContactList::remove(const int i)
{
    //cDynPrintf("clearing\n");
    // only to be used to remove persistent contacts before new contacts are added 
    assert (i < n_);
    cDynContact* c=list_[i];
    n_--;
    if (i != n_) 
    { 
        // if not at bottom of list place item at bottom of list in its place
        list_[i]=list_[n_];
        list_[i]->index_=i;
    }
    
    // remove pair if it exists
    if (c->type_ == CDYN_CONTACT && c->pair_ != NULL)
    {
        c->pair_->pair_=NULL;
        c->pair_->list_->remove(c->pair_->index_);
    }
    
    // delete any contact points associated with the contact
    if (c->type_ == CDYN_CONTACT && c->prim_ != NULL) 
    {
        while (c->head_ != NULL) 
        {
            if (!c->prim_->nodeA_->isFixed()) cDynamicsTorqueExternalZeroInwardPath(c->prim_->nodeA_);
            if (!c->prim_->nodeB_->isFixed()) cDynamicsTorqueExternalZeroInwardPath(c->prim_->nodeB_);
            delete c->head_;
        }
        // delete cDynPrimPair record
        delete(c->prim_);
        c->prim_=NULL;
    }

    delete c;
}

//---------------------------------------------------------------------------

void cDynContactList::jointLimit(const cDynTime& time, cDynJoint* joint, const cDynBoundType type)
{
    cDynContactType contact_type=(cDynContactType)type;
    // check to make sure it is not already on list
    for (int i=0;i<n_;i++) 
    {
        cDynContact* c=list_[i];
        if (c->type_ == contact_type && c->joint_ == joint)
            return;
    }
    if (time != time_)
    {
        flush();
        time_=time;
    }

    cDynContact* contact=new cDynContact;
    contact->type_=contact_type;
    contact->condition_= CDYN_UNILATERAL; // unilateral
    contact->epsilon_=joint->epsilon((cDynBoundType)contact_type);
    contact->joint_=joint;
    contact->tau_=0;
    contact->pair_=NULL;

    contact->n_=CDYN_NOTSET;

    if (n_ >= size_) grow(2*size_);
    contact->list_=this;
    contact->index_=n_;
    list_[n_++]=contact;
}
#if 0
void cDynContactList::contactUpdate(const cDynTime& time, cDynPrimPair* pair)
{
    cDynContactList* lA=pair->nodeA_->baseNode()->contact();
    cDynContactList* lB=pair->nodeB_->baseNode()->contact();

    if (time != lA->time_) 
    {
        lA->flush();
        lA->time_=time;
    }
    if (time != lB->time_)
    {
        lB->flush();
        lB->time_=time;
    }
}
#endif


//---------------------------------------------------------------------------
void cDynContactList::contact(const cDynTime& time, cDynPrimPair* pair, cDynConstraint* constraint)
{
    cDynContactList* lA=pair->nodeA_->baseNode()->contact();
    cDynContactList* lB=(pair->nodeB_->isFixed() || pair->nodeB_->baseNode()->status() == CDYN_SLEEP)?NULL:pair->nodeB_->baseNode()->contact();

    if (time != lA->time_) 
    {
        lA->flush();
        lA->time_=time;
    }
    
    if (lB != NULL && time != lB->time_) 
    {
        lB->flush();
        lB->time_=time;
    }
    
    // check to make sure it is not already on list
    for (int i=0;i<n_;i++) 
    {
        cDynContact* c=lA->list_[i];
        if (c->type_ == CDYN_CONTACT) 
        {
            if (c->prim_ == NULL) c=c->pair_;
            if (c->prim_->primA_ == pair->primA_
             && c->prim_->primB_ == pair->primB_
             && c->prim_->nodeA_ == pair->nodeA_
             && c->prim_->nodeB_ == pair->nodeB_
             && c->constraint_ == constraint) 
            {
                delete pair; // pair is not needed
                return;
            }
        }
    }
                         
    //cDynPrintf("b n_ = %d\n",n_);
    // XXXX FIX code duplicated (need proper constructors of cDynContact)
    cDynContact* cA=new cDynContact;
    cDynContact* cB=(lB != NULL)?new cDynContact:NULL;
    {
        // Create contact record for body A
        cA->type_=CDYN_CONTACT;
        cA->condition_= CDYN_UNILATERAL; // unilateral
        cA->constraint_=constraint;
        cA->prim_=pair;
        cA->numPoints_=0;
        cA->head_=NULL;
        cA->tail_=NULL;

        // add it to body A contact list
        if (lA->n_ >= lA->size_) lA->grow(2*lA->size_);
        cA->list_=lA;
        cA->index_=lA->n_;
        lA->list_[lA->n_++]=cA;
        cA->n_=CDYN_NOTSET;
        cA->pair_=cB;
    }

    if (cB != NULL)
    {
        // Create contact record for body B
        cB->type_=CDYN_CONTACT;
        cB->condition_= CDYN_UNILATERAL; // unilateral
        cB->constraint_=constraint;
        cB->prim_=NULL;		// use prim==NULL to indicate body B
        cB->numPoints_=0;
        cB->head_=NULL;
        cB->tail_=NULL;

        // add it to body B contact list
        if (lB->n_ >= lB->size_) lB->grow(2*lB->size_);
        cB->list_=lB;
        cB->index_=lB->n_;
        lB->list_[lB->n_++]=cB;
        cB->n_=CDYN_NOTSET;
        cB->pair_=cA;
    }
    //cDynPrintf("a n_ = %d\n",n_);
}

//---------------------------------------------------------------------------

void cDynContactList::setupReset()
{
    if (setup_ == CDYN_SETUP_UNINITIALIZED) return;
    setup_=CDYN_SETUP_UNINITIALIZED;
    for (int j=0;j<n_; j++) 
    {
        cDynContact* c=list_[j];
        if (c->type_ == CDYN_CONTACT && c->pair_ != NULL && c->pair_->list_->setup_ != CDYN_SETUP_UNINITIALIZED) 
        {
            c->pair_->list_->setupReset();
        }
    }
}

//---------------------------------------------------------------------------

void cDynContactList::event()
{
    if (setup_ == CDYN_SETUP_EVENTCALLBACK) return;
    setup_ = CDYN_SETUP_EVENTCALLBACK;
    for (int j=0;j<n_; j++) {
        cDynContact* c=list_[j];
        if (c->type_ == CDYN_CONTACT) 
        {
            for(cDynContactPoint* p=c->head_;p != NULL;p=p->next_) 
            {
                if (p->type_==CDYN_CONSTRAINT_NORMAL || p->type_==CDYN_CONSTRAINT_NORMAL_FRICTION) 
                {
                    if (contactEvent_ != NULL) 
                    {
                        contactEvent_->set(p,time_,1);
                        contactEvent_->handle();
                    }
                    // go through this route so that we make a callback even if contact is with fixed node
                    if (p->contact_->prim_->nodeB_->baseNode()->contact()->contactEvent_ != NULL) {
                        p->contact_->prim_->nodeB_->baseNode()->contact()->contactEvent_->set(p,time_,-1);
                        p->contact_->prim_->nodeB_->baseNode()->contact()->contactEvent_->handle();
                    }
                }
            }
        }
        if (c->type_ == CDYN_CONTACT && c->pair_ != NULL && c->pair_->list_->setup_ != CDYN_SETUP_EVENTCALLBACK)
        {
            c->pair_->list_->event();
        }
    }
}


//---------------------------------------------------------------------------

void cDynContactList::update()
{
    if (setup_ == CDYN_SETUP_UPDATE) return;
#ifdef C_DYN_DEBUG_OPS
    cDynDebugEnergy(base_,1);
#endif
    base_->discontinuity();
    base_->integrate(time_,inc_); // integrate next step
    setup_ = CDYN_SETUP_UPDATE;
    for (int j=0;j<n_; j++) 
    {
        cDynContact* c=list_[j];
        if (c->type_ == CDYN_CONTACT && c->pair_ != NULL && c->pair_->list_->setup_ != CDYN_SETUP_UPDATE) 
        {
            c->pair_->list_->update();
        }
    }
}

//---------------------------------------------------------------------------

#if 0
void cDynContactList::expand()
{
    int i;
    cDynContact* c,s;
    if (setup_ == CDYN_SETUP_EXPAND) return;

    // find points in patches
    setup_=CDYN_SETUP_EXPAND;
    for (i=0;i<n_; i++) 
    {
        c=list_[i];
        if (c->prim_) 
        {
            cA->findPoints(time_);
            if (ca->numPoints_)
                c->pair_->list_->expand();
        }
    }
    for (i=0;i<n_; i++) 
    {
        c=list_[i];
        if (c->prim_) 
        {
            for (j=0;j<i;j++)  
            {
                s=list_[i];
                if (c->pair_->list_ == s->pair_->list
            }
        }
    }
}
#endif

//---------------------------------------------------------------------------

int cDynContactList::setup(const int n)
{
    int i=n;
    bool config=false;

    if (setup_ == CDYN_SETUP_INITIALIZED) return(n);
    //base_->backup(time_);
    inc_=base_->updateInc();
#ifdef C_DYN_DEBUG_OPS
    cDynDebugEnergy(base_,0);
#endif
    base_->fwdDynamics(time_,1);

    setup_=CDYN_SETUP_INITIALIZED;
    for (int j=0;j<n_; j++) {
        cDynContact* c=list_[j];
        if (c->n_ == CDYN_NOTSET) 
        {
            c->n_=i;
            if (c->type_ == CDYN_CONTACT) 
            {
                if (c->pair_ != NULL)
                    c->pair_->n_=c->n_;

                cDynContact* cA=(c->prim_ != NULL)?c:c->pair_;
                // find contact points
                cA->findPoints(time_);
                i+= cA->numPoints_;

                if (c->pair_ != NULL) 
                {
                    if (cA->numPoints_ == 0)
                        c->pair_->n_=c->n_=CDYN_NOTSET;
                    else
                        i=c->pair_->list_->setup(i);
                } 
                else if (cA->numPoints_ == 0) c->n_=CDYN_NOTSET;
            } 
            else 
            { 
                // joint limit (1 constraint)
                i++;
            }
        }
        if (c->type_ == CDYN_CONTACT && !config) 
        {
            cDynamicsGlobalJacobian(base_->root());
            config=true;
        }
    }
    return(i);
}

//---------------------------------------------------------------------------

void cDynContactList::list(cDynContact** a)
{
    if (setup_ == CDYN_SETUP_LISTED) return;
    setup_ = CDYN_SETUP_LISTED;
    for (int i=0;i<n_; i++) 
    {
        cDynContact* c=list_[i];
        if (c->n_ != CDYN_NOTSET && a[c->n_] == NULL) 
        {
            if (c->type_ == CDYN_CONTACT) 
            {
                a[c->n_]=(c->prim_ != NULL)?c:c->pair_;
                if (c->pair_ != NULL && c->pair_->list_->setup_ != CDYN_SETUP_LISTED)
                    c->pair_->list_->list(a);
            } 
            else 
            {
                a[c->n_]=c;
                // clear status to allow true accelerations to be computed
                // XXXX CHECK IF WE STILL USE THIS c->joint_->status(CDYN_NORMAL);
            }
        }
    }
}

//---------------------------------------------------------------------------

double cDynCallbackTorqueReset(cDynJoint* joint)
{
    joint->torqueExternal(joint->jcol().dot(currentS));
    return(0.0f);
}

//---------------------------------------------------------------------------

const cDynVector3& cDynCallbackTorqueSphereReset(cDynJoint* joint)
{
    cDynVector3 tmpV0,tmpV1;
    tmpV0.multiply(joint->jsphere()[0], currentS[0]);
    tmpV1.multiply(joint->jsphere()[1], currentS[1]);
    tmpV0 += tmpV1;
    joint->storqueExternal(tmpV0);
    return(cVector3Zero);
}

//---------------------------------------------------------------------------

void cDynContactList::updatePersistent(const cDynTime time)
{
    setup_=CDYN_SETUP_UNINITIALIZED;
    time_=time;
    /*
    if (n_)
        cDynPrintf("updatePersistent (%d)\n",n_);
    else
        cDynPrintf("*(%d)\n",n_);
        */
    while (n_ > 0) 
    {
        cDynContact* c=list_[0];
        if (c->type_ == CDYN_CONTACT) 
        {
            //time_=c->prim_->checkPrimTime( , time_);
            // currently we remove all of these contacts and find them again when we need them
            remove(0);
        } 
        else 
        { 
            // joint limit
            // XXX NEED TO CHECK IF WE NEED THIS ANYMORE c->joint_->status(CDYN_NORMAL);
            c->joint_->torqueExternalZero();
            remove(0);
        }
    }
}

//---------------------------------------------------------------------------

double cDynCallbackContactIn(cDynJoint* joint)
{
    int i=currentNum;
    cDynContactPoint* p=currentContactPoint;
    double v=joint->v();
    double x= (p->s_.dot(joint->jcol()))*currentContactPointDirection;
    currentX[i] += x*v;
    return(x);
}

//---------------------------------------------------------------------------

const cDynVector3&  cDynCallbackContactSphereIn(cDynJoint* joint)
{
    static cDynVector3 x;
    int i=currentNum;
    cDynContactPoint* p=currentContactPoint;
    cDynVector3 v;
    joint->rotateHome2Local(v, joint->sv());

    cDynVector6 t;
    t[0].transposedMultiply(joint->jsphere()[0], p->s_[0]);
    t[1].transposedMultiply(joint->jsphere()[1], p->s_[1]);
    x.add(t[0], t[1]);
    x*=currentContactPointDirection;
    currentX[i] += x.dot(v);
    return(x);
}

//---------------------------------------------------------------------------

/*
const double cDynCallbackContactIn(cDynJoint* joint)
{
    int i=currentNum;
    cDynContact* c=currentContact[i];
    double v=joint->v();
    if ((double)fabs(v) > CDYN_ERROR_BOUND) {
        for (cDynContactPoint* p=c->points_;p != NULL; p=p->next_) {
            currentX[i]=p->s_.dot(joint->jcol())*v;
        }
    }
    double x= p->s_.dot(joint->jcol());
    currentX[i] += x*v;
    return(x);
}
*/

//---------------------------------------------------------------------------

double cDynCallbackLimitIn(cDynJoint* joint)
{
    cDynContact* c=currentContact[currentNum];
    assert(c->type_ != CDYN_CONTACT);
    if (c->joint_ == joint)
        return((c->type_ == CDYN_JOINT_LOWER)?1.0f:-1.0f);
    else
        return(0.0f);
}

//---------------------------------------------------------------------------

const cDynVector3&  cDynCallbackcDynmmySphereIn(cDynJoint* joint)
{
    return(cVector3Zero);
}

//---------------------------------------------------------------------------

void cDynCallbackLambdaOut(cDynJoint* joint,const double x)
{
    joint->tmp()=x;
}

//---------------------------------------------------------------------------

void cDynCallbackLambdaSphereOut(cDynJoint* joint,const cDynVector3& x)
{
    joint->tmp3()=x;
}

//---------------------------------------------------------------------------

double cDynCallbackLambdaUpdate(cDynJoint* joint)
{
    int i=currentLambdaNum;
    cDynContact* c=currentContact[i];
    assert (c->type_ == CDYN_CONTACT);
    double x=joint->tmp()*currentContactPointDirection;
    for (cDynContactPoint* p=c->head_; p != NULL && i <= currentNum; p=p->next_) 
    {
        lcol(i) += p->s_.dot(joint->jcol())*x;
        i++;
    }
    return(0.0f);
}

//---------------------------------------------------------------------------

const cDynVector3& cDynCallbackLambdaSphereUpdate(cDynJoint* joint)
{
    int i=currentLambdaNum;
    cDynContact* c=currentContact[i];
    assert (c->type_ == CDYN_CONTACT);
    cDynVector3 tmpV0,tmpV1;
    for (cDynContactPoint* p=c->head_; p != NULL && i <= currentNum; p=p->next_) 
    {
        tmpV0.multiply(joint->jsphere()[0], joint->tmp3());
        tmpV1.multiply(joint->jsphere()[1], joint->tmp3());
        lcol(i) += (p->s_[0].dot(tmpV0) + p->s_[1].dot(tmpV1))*currentContactPointDirection;
        //lcol(i) += joint->tmp3().dot(joint->jsphere()[0]*p->s_[0] 
        //	 +  joint->jsphere()[1]*p->s_[1])*currentContactPointDirection;
        i++;
    }
    return(cVector3Zero);
}

//---------------------------------------------------------------------------

void cDynContactList::updateLambda(double direction)
{
    currentContactPointDirection = direction;
    cDynContact* current=NULL;
    if (direction == 0) current=currentContact[currentNum];
    else current=currentContactPoint->contact_;
    cDynBaseNode* base=NULL;
    cDynObject* obj=NULL;
    if (current->type_ == CDYN_CONTACT) 
    {
        if (direction < 0) 
        {
            base=current->list_->base_;
            obj=current->prim_->nodeA_;
        } 
        else 
        {
            base=(current->pair_ == NULL)?NULL:current->pair_->list_->base_;
            obj=(current->pair_ == NULL)?NULL:current->prim_->nodeB_;
        }
        if (base == NULL) return; // contact with fixed geometry
//		base->fwdCallbacks(&cDynCallbackContactIn, &DeCallbackContactSphereIn,
//			&cDynCallbackLambdaOut, &cDynCallbackLambdaSphereOut);
        cDynamicsFwdDynamicsConfig(base->root(), obj, &cDynCallbackContactIn, &cDynCallbackContactSphereIn,
            &cDynCallbackLambdaOut, &cDynCallbackLambdaSphereOut, 0);

    } 
    else 
    {
        base=current->list_->base_;
        obj=current->joint_->object();
//		base->fwdCallbacks(
//			&DeCallbackLimitIn, NULL,
//			&cDynCallbackLambdaOut, &cDynCallbackLambdaSphereOut);
        cDynamicsFwdDynamicsConfig(base->root(), obj, &cDynCallbackLimitIn, &cDynCallbackcDynmmySphereIn,
            &cDynCallbackLambdaOut, &cDynCallbackLambdaSphereOut, 0);
    }
    
//	base->fwdCallbacks(&DeCallbackLambdaUpdate, &DeCallbackLambdaSphereUpdate, NULL, NULL);
    int i=0;
    while (i <= currentNum) 
    {
        cDynContact* c=currentContact[i];
        if (c->type_ == CDYN_CONTACT) {
            if (!c->prim_->nodeA_->isFixed() && base == c->prim_->nodeA_->baseNode()) 
            {
                currentLambdaNum=i;
                currentContactPointDirection = -1.0f;
                cDynamicsFwdDynamicsIn(c->prim_->nodeA_,&cDynCallbackLambdaUpdate, &cDynCallbackLambdaSphereUpdate);
            }
            if (!c->prim_->nodeB_->isFixed() && base == c->prim_->nodeB_->baseNode()) 
            {
                currentLambdaNum=i;
                currentContactPointDirection = 1.0f;
                cDynamicsFwdDynamicsIn(c->prim_->nodeB_,&cDynCallbackLambdaUpdate, &cDynCallbackLambdaSphereUpdate);
            }
            i += c->numPoints_;
        } 
        else 
        {
            if (base == c->list_->base_) 
            {
                if (c->type_ == CDYN_JOINT_LOWER)
                    lcol(i) += c->joint_->tmp();
                else 
                    lcol(i) -= c->joint_->tmp();
            }
            i++;
        }
    }
}
#if 0

static const double CDYN_LAMBDA_EPSILON = CDYN_ERROR_BOUND;

void cDynCallbackLambdaOut(cDynJoint* joint,const double x)
{
    if (fabs(x) < CDYN_LAMBDA_EPSILON)
      {
        return;
      }	
    int i=0;
    while (i <= currentNum) 
    {
        cDynContact* c=currentContact[i];
        if (c->type_ == CDYN_CONTACT) 
        {
            double dx=x*currentContactPointDirection;
            cDynContactPoint* p=c->head_;
            while (p != NULL && i <= currentNum) 
            {
                lcol(i) += p->s_.dot(joint->jcol())*dx;
                p=p->next_; i++;
            }
        } 
        else 
        {
            if (joint == c->joint_) 
            {
                if (c->type_ == CDYN_JOINT_LOWER)
                    lcol(i) += x;
                else
                    lcol(i) -= x;
            }
            i++;
        }
    }
}

void cDynCallbackLambdaSphereOut(cDynJoint* joint,const cDynVector3& x)
{
    int i=0;
    while (i <= currentNum) 
    {
        cDynContact* c=currentContact[i];
        if (c->type_ == CDYN_CONTACT) 
        {
            for (cDynContactPoint* p=c->head_; p != NULL; p=p->next_) 
            {
                lcol(i) += (p->s_[0].dot(joint->jsphere()[0]*x)
                  + p->s_[1].dot(joint->jsphere()[1]*x))*currentContactPointDirection;
                i++;
            }
        } 
        else 
        {
            i++;
        }
    }
}
#endif

//---------------------------------------------------------------------------

double cDynCallbackJtIn(cDynJoint* joint)
{
    int i=currentNum;
    cDynContact* c=currentContact[i];
    double x;
    if (c->type_ == CDYN_CONTACT) 
    {
        x=0.0f;
        for (cDynContactPoint* p=c->head_; p != NULL; p=p->next_) 
        {
            x +=p->s_.dot(joint->jcol())*currentX[i]*currentContactPointDirection;
            i++;
        }
    } 
    else 
    {
        if (joint == c->joint_) 
        {
            x=(c->type_ == CDYN_JOINT_LOWER)?currentX[i]:-currentX[i];
        } 
        else 
        {
            x=0.0f;
        }
    }
    
    if (pass == CDYN_FORCE) 
    {
        joint->torqueExternal(x);
    }
    return(x);
}

//---------------------------------------------------------------------------

const cDynVector3& cDynCallbackJtSphereIn(cDynJoint* joint)
{
    int i=currentNum;
    cDynContact* c=currentContact[i];
    double y;
    if (c->type_ == CDYN_CONTACT) 
    {
        retV.zero();
        for (cDynContactPoint* p=c->head_; p != NULL; p=p->next_) 
        {
            y=currentX[i]*currentContactPointDirection;
            if (y != 0.0f) 
            {
                cDynVector6 t;
                t[0].transposedMultiply(joint->jsphere()[0], p->s_[0]);
                t[1].transposedMultiply(joint->jsphere()[1], p->s_[1]);
                cDynVector3 tmpV;
                tmpV.add(t[0], t[1]);
                tmpV*=y;
                retV += tmpV;
            }
            i++;
        }
        if (pass == CDYN_FORCE) 
        {
            joint->storqueExternal(retV);
        }
        return(retV);
    } 
    else 
    {
        return(cVector3Zero);
    }
}

//---------------------------------------------------------------------------

void cDynCallbackPositionOut(cDynJoint* joint, const double x)
{
    joint->qerr(x);
}

//---------------------------------------------------------------------------

void cDynCallbackPositionSphereOut(cDynJoint* joint, const cDynVector3& x)
{
    joint->sqerr(x);
}

//---------------------------------------------------------------------------

void cDynCallbackVelocityOut(cDynJoint* joint, const double x)
{
    joint->v(joint->v()+x);
}

//---------------------------------------------------------------------------

void cDynCallbackVelocitySphereOut(cDynJoint* joint, const cDynVector3& x)
{
    cDynVector3 r = x;
    r += joint->sv();
    joint->sv(r);
}

//---------------------------------------------------------------------------

void cDynCallbackAccelerationOut(cDynJoint* joint, const double x)
{
    joint->a(joint->a()+x);
}

//---------------------------------------------------------------------------

void cDynCallbackAccelerationSphereOut(cDynJoint* joint, const cDynVector3& x)
{
    cDynVector3 r = x;
    r += joint->sa();
    joint->sa(r);
}

//---------------------------------------------------------------------------

double cDynCallbackAfreeIn(cDynJoint* joint)
{
    int i=currentNum;
    cDynContact* c=currentContact[i];
    int ii=c->n_;
    assert (c->type_ == CDYN_CONTACT);
    //cDynVector6 v;
    cDynVector6 tmpV;
    for (cDynContactPoint* p=c->head_; p != NULL; p=p->next_) 
    {
        //c->list_->base_->dynamics_->globalBiasAcceleration(c->node_,&v);
        //currentX[ii] +=p->s_.dot(joint->jcol()*joint->a() + v);
        tmpV.multiply(joint->jcol(), joint->a());
        currentX[ii] += p->s_.dot(tmpV)*currentContactPointDirection;
        i++;ii++;
    }
    return(0);
}

//---------------------------------------------------------------------------

const cDynVector3& cDynCallbackAfreeSphereIn(cDynJoint* joint)
{
    int i=currentNum;
    cDynContact* c=currentContact[i];
    assert (c->type_ == CDYN_CONTACT);
    int ii=c->n_;
    //cDynVector6 v;
    cDynVector3 tmpV;
    for (cDynContactPoint* p=c->head_; p != NULL; p=p->next_)
    {
        //c->list_->base_->dynamics_->globalBiasAcceleration(c->node_,&v);
        //currentX[ii] += p->s_[0].dot(joint->jsphere()[0]*joint->sa()+v[0]);
        //currentX[ii] += p->s_[1].dot(joint->jsphere()[1]*joint->sa()+v[1]);
        tmpV.multiply(joint->jsphere()[0], joint->sa());
        currentX[ii] += p->s_[0].dot(tmpV)*currentContactPointDirection;
        tmpV.multiply(joint->jsphere()[1], joint->sa());
        currentX[ii] += p->s_[1].dot(tmpV)*currentContactPointDirection;
        i++;ii++;
    }
    return(cVector3Zero);
}

//---------------------------------------------------------------------------

void cDynContactList::force(const int m, cDynContact** contact, double* f)
{
    for (int i=0;i<m;) 
    {
        cDynContact* c=contact[i];
        int ii=c->n_;
        if (c->type_ == CDYN_CONTACT) 
        {
            for (cDynContactPoint* p=c->head_;p != NULL; p=p->next_) 
            {
                f[ii++] += p->force_;i++;
            }
        } 
        else 
        {
            f[ii] += c->tau_;i++;
        }
    }
}

//---------------------------------------------------------------------------

void cDynContactList::acceleration(const int m, cDynContact** contact, double* a) 
{
    currentX=a;
    for (int i=0;i<m;) 
    {
        cDynContact* c=contact[i];
        int ii=c->n_;
        if (c->type_ == CDYN_CONTACT) 
        {
            for (int j=0;j<c->numPoints_;j++)
                a[ii+j]=0.0f;
            cDynContact* p=c->pair_;
            currentNum=i;
//			c->list_->base_->fwdCallbacks(&DeCallbackAfreeIn, &DeCallbackAfreeSphereIn, NULL, NULL);
            if (c->list_->setup_ != CDYN_SETUP_ACCELERATION) 
            {
                cDynamicsFwdDynamicsVel(c->list_->base_->root());
//				cDynamicsFwdDynamics(c->list_->base_->root(),&c->list_->base_->gravity());
                c->list_->setup_ = CDYN_SETUP_ACCELERATION;
            }
            currentContactPointDirection = -1.0f;
            cDynamicsFwdDynamicsIn(c->prim_->nodeA_, &cDynCallbackAfreeIn, &cDynCallbackAfreeSphereIn);
            if (p != NULL) 
            {
//				p->list_->base_->fwdCallbacks(&DeCallbackAfreeIn, &DeCallbackAfreeSphereIn, NULL, NULL);
                if (p->list_->setup_ != CDYN_SETUP_ACCELERATION) 
                {
                    cDynamicsFwdDynamicsVel(p->list_->base_->root());
//					cDynamicsFwdDynamics(p->list_->base_->root(),&p->list_->base_->gravity());
                    p->list_->setup_ = CDYN_SETUP_ACCELERATION;
                }
                currentContactPointDirection =  1.0f;
                cDynamicsFwdDynamicsIn(c->prim_->nodeB_,&cDynCallbackAfreeIn, &cDynCallbackAfreeSphereIn);
            }
            i += c->numPoints_;
        } 
        else 
        { 
            // joint limit
            if (c->list_->setup_ != CDYN_SETUP_ACCELERATION) 
            {
                cDynamicsFwdDynamicsVel(c->list_->base_->root());
//				cDynamicsFwdDynamics(c->list_->base_->root(),&c->list_->base_->gravity());
                c->list_->setup_ = CDYN_SETUP_ACCELERATION;
            }
            
            if (c->type_ == CDYN_JOINT_LOWER)
                a[ii] = c->joint_->a();
            else
                a[ii] = -c->joint_->a();
            i++;
        }
    }
}

//---------------------------------------------------------------------------
#define NCOL 10
//---------------------------------------------------------------------------

void cDynDebugResolve(int m, int n,int index[], int type[][3],double b[], double x[])
{
    // print stuff
    for (int j=0;j<m;j+=NCOL) 
    {
    cDynPrintf("set   ");
    {for (int i=j;i<m && i<j+NCOL;i++) 
    {
        char c=(i<n)?'C':'_';
        char t='E';
        switch (BASE(type[index[i]][pass])) 
        {
            case CDYN_UNILATERAL: t= ' '; break;
            case CDYN_BILATERAL:  t= '*'; break;
            case CDYN_REDUNDANT:  t= 'R'; break;
            case CDYN_NOTACTIVE:  t= '#'; break;
        }
        char f=' ';
        switch (FLAG(type[index[i]][pass])) 
        {
            case CDYN_FRICTION_NORMAL: f= '^'; break;
            case CDYN_SATURATED_NORMAL: f= 'S'; break;
            case CDYN_FRICTION: f= 'F'; break;
        }
        cDynPrintf("%10d%c%c%c ",index[i], c, f, t);
    }}
    
    cDynPrintf("\nb     ");
    {for (int i=j;i<m && i<j+NCOL;i++) 
        cDynPrintf("%13.9f ", b[index[i]]);
    }
    
    cDynPrintf("\nx     ");
    {for (int i=j;i<m && i<j+NCOL;i++) 
        cDynPrintf("%13.9f ", x[index[i]]);
    }

    cDynPrintf("\n");
    }
}

//---------------------------------------------------------------------------

bool cDynContactList::resolve(const int m)
{

    int *index = (int *)CDYN_ALLOCA(m*sizeof(int));
    int *rindex = (int*)CDYN_ALLOCA(m*sizeof(int));

    double* Lambda = NULL;
    double* diag = NULL;
    int l=m*m;
    if ((m*((m+5)*sizeof(double)+2*sizeof(cDynContact *))) < MAXSTACKALLOC) 
    {
        Lambda=(double *)CDYN_ALLOCA(m*((m+5)*sizeof(double)+2*sizeof(cDynContact *)));
        diag = Lambda+l;
    }

    bool bigLambda=false;
    if (Lambda == NULL) 
    {
        Lambda= (double *)malloc(l*sizeof(double));
        if (cDynDebugLevel == 0)
            cDynPrintf("cDynContactList:resolve: malloc: %d double\n",l);
        diag = (double *)CDYN_ALLOCA(m*(5*sizeof(double)+2*sizeof(cDynContact *)));
        bigLambda=true;
    }

    double* b = diag+m;
    double* x = b+m;
    double* err = x+m;
    double* maxerr = err+m;
    cDynContact** contact = (cDynContact **)(maxerr+m);
    cDynContactPoint** contactPoints = (cDynContactPoint**)(contact+m);

    double* a= (double *)CDYN_ALLOCA(m*sizeof(double));
    double* impulse= (double *)CDYN_ALLOCA(m*sizeof(double));
    double* force= (double *)CDYN_ALLOCA(m*sizeof(double));
    double* e= (double *)CDYN_ALLOCA(m*sizeof(double));
    double* inc= (double *)CDYN_ALLOCA(m*sizeof(double));
    double* v= (double *)CDYN_ALLOCA(m*sizeof(double));
    double* dv= (double *)CDYN_ALLOCA(m*sizeof(double));
    double* delta= (double *)CDYN_ALLOCA(m*sizeof(double));
    double* max = (double*)CDYN_ALLOCA(m * sizeof(double));
    int (*type)[3] = (int (*)[3])CDYN_ALLOCA(3*m*sizeof(int));

    //cDynPrintf("setting start(%d)\n",n_);
    testB=b;
    testE=e;
    testErr=err;
    testX=x;
    testV=v;
    testDV=dv;
    testA=a;
    testDiag=diag;
    testDelta=delta;
    testMax=max;
    //testType=type;
    static long iter=0;
    static long rowiter=0;
    static int flag=0;
    int cImpulse=0;
    int cForce=0;
    int cError=0;
    int cZero=0;
    int cN=0;
    int cNV=0;

    flag=(cDynDebugLevel > 2)?1:0;
    if (cDynDebugLevel > 1)
    cDynPrintf("iter=%ld t=%5.9f m=%d\n",iter,time_,m);
    iter++;

    if (bigLambda)
    {
        memset(Lambda,0,l*sizeof(double));
        memset(diag,0,m*(5*sizeof(double)+2*sizeof(cDynContact *)));
    }
    else
    {
        memset(Lambda,0,m*((m+5)*sizeof(double)+2*sizeof(cDynContact *)));
    }

    // build up a list of contacts in an array for easy access
    list(contact);

    //-----------------------------------------------------------------------
    // DEBUG STUFF START
    //-----------------------------------------------------------------------

    if (flag) 
    {
        {for (int i=0;i<m;) 
        {
            cDynContact *c=contact[i];
            if (c->type_ == CDYN_CONTACT) 
            {
                cDynPrintf("%s(%d) <---> %s(%d)\n",
                    (c->prim_->nodeA_ != NULL)?cDynObjectName(c->prim_->nodeA_):"X",
                    (c->prim_->primA_ != NULL)?c->prim_->primA_->id:0,
                    (c->prim_->nodeB_ != NULL)?cDynObjectName(c->prim_->nodeB_):"X",
                    (c->prim_->primB_ != NULL)?c->prim_->primB_->id:0
                    );
                for(cDynContactPoint* p=c->head_;p != NULL;p=p->next_) 
                {
                    char c=' ';
                    switch (p->type_) 
                    {
                    case CDYN_CONSTRAINT_NORMAL: c='n'; break;
                    case CDYN_CONSTRAINT_NORMAL_FRICTION: c='N'; break;
                    case CDYN_CONSTRAINT_FRICTION_X: c='x';break;
                    case CDYN_CONSTRAINT_FRICTION_Y: c='y';break;
                    default: c='*';
                    }
                    cDynPrintf("\t%d: %c=[%7.5f,%7.5f,%7.5f] r=[%7.5f,%7.5f,%7.5f] err=%5.5e\n",i,
                        c,
                        p->s_[0][0], p->s_[0][1], p->s_[0][2],
                        p->r_[0], p->r_[1], p->r_[2],
                        p->err_
                        );

                    i++;
                }
            } 
            else 
            {
                double err;
                char ctype;
                if (c->type_ == CDYN_JOINT_UPPER)
                {
                    err=c->joint_->bound(CDYN_UPPER) - c->joint_->q();
                    ctype='U';
                } 
                else 
                {
                    err=c->joint_->q() - c->joint_->bound(CDYN_LOWER);
                    ctype='L';
                }

                cDynPrintf("\t%d: joint (%d) %c %s (%s) err=%5.5e\n",i,c->joint_->object()->uid(),ctype,c->joint_->data(),cDynObjectName(c->joint_->object()),err);
                i++;
            }
        }}
        cDynPrintf("\n");
    }

    //-----------------------------------------------------------------------
    // DEBUG END
    //-----------------------------------------------------------------------

    //make links to some data structures available globally for
    //the callback functions
    currentSize=m;
    currentContact=contact;
    currentLambda=Lambda;

    testIndex=index;
    {for (int i=0;i<m;i++) index[i]=rindex[i]=i;}
    currentIndex=index;

#ifdef RTIMERLAMBDA
    rtimer.start();
#endif
    
    //-----------------------------------------------------------------------
    //
    // INITIALIZE LAMBDA, VELOCITY(v)
    //
    //-----------------------------------------------------------------------

    {for (int i=0; i<m; ) 
    {
        cDynContact* c=contact[i];

        //----------------------------------------------------------------------- 
        // CDYN_CONTACT
        //-----------------------------------------------------------------------

        if (c->type_ == CDYN_CONTACT) 
        {
            double tinc=(c->pair_ == NULL)?
            c->list_->inc_:
            ((c->list_->inc_ < c->pair_->list_->inc_)?c->list_->inc_:c->pair_->list_->inc_);
            currentX=v;
    
            // compute row of Lambda for each contact point
            for(cDynContactPoint* p=c->head_;p != NULL;p=p->next_)
            {
                int itype=0,ftype=0,etype=0;
                switch (p->type_) 
                {
                    case CDYN_CONSTRAINT_NORMAL: itype = ftype = etype = CDYN_UNILATERAL; break;
                    case CDYN_CONSTRAINT_BILATERAL: itype = ftype = etype = CDYN_BILATERAL; break;
                    case CDYN_CONSTRAINT_NORMAL_FRICTION: itype = ftype = etype = CDYN_FRICTION_NORMAL + CDYN_UNILATERAL; break;
                    case CDYN_CONSTRAINT_FRICTION_X:
                    case CDYN_CONSTRAINT_FRICTION_Y: itype = ftype = etype = CDYN_FRICTION + CDYN_NOTACTIVE; break;
                    default: break;
                }

                type[i][CDYN_IMPULSE]=itype;
                type[i][CDYN_FORCE]=ftype;
                type[i][CDYN_ERROR]=etype;
                contactPoints[i]=p;
                if (BASE(itype) == CDYN_UNILATERAL)
                e[i]=c->prim_->epsilon();
                else
                e[i]=0.0f;

                err[i]=p->err_;
                maxerr[i]=p->maxerr_;
                inc[i]=tinc;
                if (p->info_ != NULL) 
                {
                    v[i] = p->info_->velocity();
                    if (p->err_ == 0.0) type[i][CDYN_ERROR]=CDYN_NOTACTIVE;
                }
                else 
                    v[i]=0.0f;

        #undef CDYN_FRICTION_WAIT
        #ifdef CDYN_FRICTION_WAIT
                if (itype != CDYN_FRICTION + CDYN_NOTACTIVE)
        #endif
                { 
                    // dont compute lambda for friction constraints until needed below
                    currentContactPoint=p;
                    currentNum=i;
                    currentCol=m*i;

                    updateLambda(-1.0f);
                    updateLambda( 1.0f);

                    if (lcol(i) < 1e-12) 
                    { 
                        // degenerate contact, don't need to consider from now on
                        //base_->world()->ignore(c->prim_->nodeA_,c->prim_->nodeB_);
                    }

        #ifdef C_DYN_DEBUG_OPS
                    vminus[i]=v[i];
        #endif

                    if (BASE(type[i][CDYN_IMPULSE]) == CDYN_BILATERAL || (v[i] > 0.0f && -v[i]*e[i] < VLIMIT)) 
                    {
                        b[i]=v[i];
                        v[i]=0.0f;
                    } 
                    else 
                    {
                        b[i]= (1.0f+e[i])*v[i];
                        v[i] *= -e[i];
                    }
                }
                i++;
            }
        } 

        //----------------------------------------------------------------------- 
        // CDYN_JOINT_LOWER / CDYN_JOINT_UPPER
        //-----------------------------------------------------------------------

        else 
        { 
            // joint limit
            // determine the delta T
            {for (int j=CDYN_IMPULSE;j<CDYN_DONE;j++)
                type[i][j]=c->condition_;
            }
            
            if (BASE(type[i][CDYN_IMPULSE]) == CDYN_UNILATERAL)
                e[i]=c->epsilon_;
            else
                e[i]=0.0f;

            if (c->type_ == CDYN_JOINT_UPPER) 
            {
                err[i]=(c->joint_->bound(CDYN_UPPER) - c->joint_->q());
                maxerr[i]=c->joint_->error(CDYN_UPPER);
            } 
            else 
            {
                err[i]=(c->joint_->q() - c->joint_->bound(CDYN_LOWER));
                maxerr[i]=c->joint_->error(CDYN_LOWER);
            }

            inc[i]=c->list_->inc_;

            currentNum=i;
            currentCol=m*i;
            updateLambda(0);
      
            double vel=c->joint_->v();
            if (c->type_ == CDYN_JOINT_UPPER) vel = -vel;
#ifdef C_DYN_DEBUG_OPS
            vminus[i]=vel;
#endif
            v[i]= -e[i]*vel;

            if (e[i] == 0.0f || (vel > 0.0f && v[i] < VLIMIT)) 
            { 
                // disappearing contact
                b[i]=vel;
                v[i]=0.0;
            } 
            else 
            {
                b[i]= (1.0f+e[i])*vel;
            }

            assert(v[i]>=0.0f);

            i++;
        }
    }
    }

    // test adding small value to lambda
    // HOO
    for (int i=0;i<m;i++) 
    {
        lambda(i,i) += 1e-6;
    }

    //-----------------------------------------------------------------------
    //
    // END INITIALIZATION
    //
    //-----------------------------------------------------------------------

#ifdef RTIMERLAMBDA
    rtimer.end();
    rtotal += rtimer.time();
    rcount++;
#endif


//---------------------------------------------------------------------------
#ifdef CDYN_DEBUG  // BEGIN DEBUG STUFF
//---------------------------------------------------------------------------
    
    if (flag) 
    {
        display(m,Lambda,b);
        displaySmall(m,Lambda,b);

    //-----------------------------------------------------------------------
    #if 0
    //-----------------------------------------------------------------------
        {for (int i=0;i<m;i++) 
        {
            if (contact[i]->type_ != CDYN_CONTACT) 
            {
                cDynPrintf("%s [%12.8f,%12.8f,%12.8f]\n",
                    contact[i]->joint_->data(),
                    contact[i]->joint_->q(),
                    contact[i]->joint_->v(),
                    contact[i]->joint_->a()
                    );
            } 
            else 
            {
                cDynPrintf("%s(%d) <---> (%d)%s\n",
                    cDynObjectName(contact[i]->prim_->nodeA_),
                    contact[i]->prim_->nodeA_->id(),
                    contact[i]->prim_->nodeB_->id(),
                    cDynObjectName(contact[i]->prim_->nodeB_)
                    );
            }
        }
        }
        
        if (m > 0) cDynPrintf("\n");
        base_->root()->displaySmallTree(NULL);

    //-----------------------------------------------------------------------
    #endif
    //-----------------------------------------------------------------------
    }
   
//---------------------------------------------------------------------------
#endif // END DEBUG STUFF
//---------------------------------------------------------------------------

    //-----------------------------------------------------------------------
    //
    // BEGIN SOLUTION OF LCP
    //
    //-----------------------------------------------------------------------

    int n=0;
    int nv=0;
    for (pass=CDYN_IMPULSE; pass < CDYN_DONE; pass++) 
    {
        n=0; nv=0;
        bool done=false;

        #ifdef CHECKMATH
        // DEBUG functions
        double debug[m];
        for (int k=0;k<m;k++) debug[k]=b[k];
        #endif

        if (flag) 
        {
            switch (pass) 
            {
            case CDYN_IMPULSE: cDynPrintf("IMPULSE PASS\n");break;
            case CDYN_FORCE: cDynPrintf("FORCE PASS\n");break;
            case CDYN_ERROR: cDynPrintf("ERROR PASS\n");break;
            }

            cDynDebugResolve(m,n,index,type,b,x);
        }

        int loopCount = 0; // make sure loopCount ~ m
        int dindex=-1;
        int addfriction=0;

        while (!done) 
        {
            if (pass == CDYN_FORCE) 
            { 
                // determine which constraints to consider for the force pass
                for (int i=n;i<m;i++) 
                {	
                    // UPDATE NOTACTIVE CONSTRAINTS FOR FORCE PASS
                    // if constraint will be free within one time cycle do not consider it during force pass
                    if (i != dindex && BASE(type[index[i]][CDYN_IMPULSE]) == CDYN_UNILATERAL) 
                    {
                        type[index[i]][CDYN_FORCE]= FLAG(type[index[i]][CDYN_FORCE]) + CDYN_UNILATERAL;
            
                        if (v[index[i]] > 0.0f || b[index[i]] > 0.0f) 
                        {
                            double h=(err[index[i]] > 0.0f)?maxerr[index[i]]-err[index[i]]:maxerr[index[i]];
                            double sq=v[index[i]]*v[index[i]] + 2.0f*b[index[i]]*h;
                            if (sq > 0.0f) 
                            {
                                if (v[index[i]] >= 0.0f) 
                                {
                                    double q= -0.5f*(v[index[i]]+cDynSqrt(sq)); // q must be leq 0
                                    if ((b[index[i]] < -1e-8f && q/b[index[i]] < inc[index[i]])
                                         || (q < -1e-8f && -h/q <inc[index[i]]))
                                    {
                                        type[index[i]][CDYN_FORCE]=FLAG(type[index[i]][CDYN_FORCE]) + CDYN_NOTACTIVE;
                                    }
                                } 
                                else 
                                {
                                    double q= -0.5f*(v[index[i]]-cDynSqrt(sq)); // q must be geq 0
                                    if (b[index[i]] > 1e-8f && q/b[index[i]] < inc[index[i]])
                                    {
                                        type[index[i]][CDYN_FORCE]=FLAG(type[index[i]][CDYN_FORCE]) + CDYN_NOTACTIVE;
                                    }
                                }
                            }
                            
                            if (flag && BASE(type[index[i]][CDYN_FORCE])==CDYN_NOTACTIVE) 
                                cDynPrintf("constraint %d force tagged as NOT ACTIVE\n",index[i]);
                        }
                    }
                }
            }

            // Find next constraint to consider
            double min= -1.0e-8f;
            //double min= -1.0e-6f;
            int mindex=-1;
        
            if (dindex != -1 && (BASE(type[index[dindex]][pass]) == CDYN_BILATERAL || b[index[dindex]] < 0.0f)) 
            { 
                // variable is already being driven to zero, continue
                mindex=dindex;
            } 
            else 
            {
                for (int i=n;i<m;i++) 
                {
                    if (BASE(type[index[i]][pass]) == CDYN_BILATERAL) 
                    {
                        mindex=i;
                        if (contactPoints[index[i]]->info_ != NULL) 
                        {
                            switch (pass) 
                            {
                                case CDYN_IMPULSE: max[index[i]]=contactPoints[index[i]]->info_->maxImpulse(); break;
                                case CDYN_FORCE: max[index[i]]=contactPoints[index[i]]->info_->maxForce(); break;
                                default: max[index[i]]=0.0;
                            }
                        } 
                        else if (type[index[i]][pass] == CDYN_BILATERAL) 
                        { 
                            // user defined constraint
                            max[index[i]]=0.0;
                        }
                        break;
                    } 
                    else if (BASE(type[index[i]][pass]) == CDYN_UNILATERAL && b[index[i]] < min) 
                    {
                        min=b[index[i]];
                        mindex=i;
                    }
                }
            }
      
            if (mindex == -1 && !addfriction && pass != CDYN_ERROR) 
            {
                addfriction++;
                for (int i=n;i<m;i++) 
                {
                    if (type[index[i]][pass] == CDYN_FRICTION + CDYN_NOTACTIVE) 
                    {
                        int k=index[i]-1;
                        if (FLAG(type[k][pass]) != CDYN_FRICTION_NORMAL) 
                            k--;
                
                        double umax=contactPoints[k]->fr_->ustatic()*x[k];
                
                        // ADD GRIP FRICTION (this is a hack because we dont have the normal force yet at the time we compute the impulses)
                
                        if (v[k] <= 0.0 && umax < contactPoints[k]->fr_->ugrip())
                        umax = contactPoints[k]->fr_->ugrip();
                        if (umax > 1e-6) 
                        {
                            max[index[i]]=umax;
                            type[index[i]][pass] = CDYN_FRICTION + CDYN_BILATERAL;

                            // HANDLE DYNAMIC and VISCOUS FRICTION
                            if (pass == CDYN_FORCE) 
                            {
                                if (contactPoints[k]->fr_->udynamic() != 0.0f) 
                                {
                                    double vmag=0.0;
                                    for (int l=k+1;FLAG(type[l][CDYN_FORCE]) == CDYN_FRICTION; l++) 
                                    {
                                        vmag += v[l]*v[l];
                                    }
                                    if (vmag > 1e-6) 
                                    {
                                        b[index[i]] += contactPoints[k]->fr_->udynamic()*lambda(k,k)*x[k]*v[index[i]]/cDynSqrt(vmag);
                                    }
                                }
                                b[index[i]] += contactPoints[k]->fr_->uviscous()*v[index[i]];
                            }
#ifdef CDYN_FRICTION_WAIT
                            // may need to compute Lambda for new constraint
                            if (lambda(i,i,index) == 0.0) 
                            {
                                int j=index[i];
                                currentContactPoint=contactPoints[j];
                                currentNum=j;
                                currentCol=m*j;

                                updateLambda(-1.0f);
                                updateLambda( 1.0f);
                                if (pass == CDYN_IMPULSE) 
                                {
                                    b[i] += v[i];
                                    v[i]=0.0f;
                                }
                            }
#endif
                            mindex=i;

                        }
                    }
                }
            }

            if (mindex == -1) 
            {
                // remaining constraints all ok with x=0
                done=true;
                {for (int i=0;i<m;i++) 
                {
                    if (BASE(type[index[i]][pass])==CDYN_NOTACTIVE)
                    {
                        assert(i>=n);
                    }
                    if (BASE(type[index[i]][pass])==CDYN_BILATERAL)
                    {
                        assert(i<n);
                    }
                }
                }
          
                switch (pass) 
                {
                    case CDYN_IMPULSE: cImpulse = n; break;
                    case CDYN_FORCE:   cForce = n; break;
                    case CDYN_ERROR:   cError = n; break;
                }
            } 
            else 
            {
                loopCount++;

                rowiter++;
                if (flag) cDynPrintf("Drive Variable %d (%ld)\n",index[mindex],rowiter);

                // swap worst offender to next spot
    
                int tmp=index[n];
                index[n]=index[mindex];
                index[mindex]=tmp;

                rindex[index[n]]=n;
                rindex[index[mindex]]=mindex;

#if 1           // original + opt
                int k=nv;
                while (nv <= n) 
                {
                    if (nv==n)
                        k=0;

                    // do a forward subsitution. This is part of an
                    // incremental cholesky solution and used to compute
                    // delta.
                    {for (int j=k;j<nv;j++) 
                    {
                        double ll=lambda(j,nv,index);
                        for (int i=0;i<j;i++)
                            ll -= L(j,i,index)*L(nv,i,index);
                        L(nv,j,index) = ll/diag[index[j]];
                    }}

                    double dd=lambda(nv,nv,index);
                    {for (int i=0;i<nv;i++) 
                    {
                        double l=L(nv,i,index);
                        dd -= l*l;
                    }}

                    // check if matrix is singular
                    // CHANGE
                    if (dd <= 1e-9f) 
                    { 
                        if (x[index[nv]] != 0.0f) 
                        { 
                            // if constraint is redundant with respect to current constraints redistribute force/impulses
                            // do a backward substition
                            for (int j=nv-1;j>=0;j--) 
                            {
                                delta[j] = L(nv,j,index);
                                for (int i=j+1;i<nv;i++)
                                    delta[j] -= L(i,j,index)*delta[i];
                                delta[j] /= diag[index[j]];
                            }
                            
                            double u= x[index[nv]];
                            for (int i=0;i<nv;i++) 
                            {
                                x[index[i]] += u*delta[i];
                            }
                            x[index[nv]]=0.0f;
                        }
            
                        if (BASE(type[index[nv]][pass]) == CDYN_BILATERAL) 
                        {
                            type[index[nv]][pass] = FLAG(type[index[nv]][pass]) + CDYN_REDUNDANT;
                            if (flag) cDynPrintf("constraint %d marked a redundant\n",index[nv]);
                        } 
                        else 
                        {
                            if (flag) cDynPrintf("constraint %d redistributed %d\n",index[nv],type[index[nv]][pass]);
                            b[index[nv]]=0.0f;
                        }
            
                        if (n != nv) 
                        {
                            assert(n != 0);
                            n--;
                            int tmp=index[n];
                            index[n]=index[nv];
                            index[nv]=tmp;
                        }
                        break;
                    } 
                    else 
                    {
                        diag[index[nv]] = cDynSqrt(dd);
                        nv++;
                    }
                }
#endif
                // if constraint was found to be redundant then start again
                if (nv > n) 
                {

                    // BEGIN DEBUG
                    if (false) 
                    { 
                        // used during debugging
                        cDynPrintf("done: n=%d\n",n);
                        for(int i=0;i<=n;i++) 
                        {
                            for(int j=0;j<=i;j++) 
                            {
                                if (i==j)
                                    cDynPrintf("%5.4f ",diag[index[i]]);
                                else
                                    cDynPrintf("%5.4f ",L(i,j,index));
                            }
                            cDynPrintf("\n");
                        }
                    }
                    // END DEBUG

                    // backward subsitution
                    {for (int j=n-1;j>=0;j--) 
                    {
                        delta[j] = L(n,j,index);
                        for (int i=j+1;i<n;i++)
                            delta[j] -= L(i,j,index)*delta[i];
                        delta[j] /= diag[index[j]];
                    }}

                    // BEGIN DEBUG
                    if (false) 
                    {
                        cDynPrintf("delta ");
                        for (int i=0;i<n;i++) 
                        {
                            cDynPrintf("%12.9f ", delta[i]);
                        }
                        cDynPrintf("\n");
                    }
                    // END DEBUG

                    //-------------------------------------------------------------------
                    //
                    // MAXSTEP SEGMENT: determine max step size
                    //
                    //-------------------------------------------------------------------

                    // first use the new point (note: we are actually computing negative acceleration)
                    delta[n]= -lambda(n,n,index);
                    {for (int i=0;i<n;i++) 
                    {
                        delta[n] += lambda(n,i,index)*delta[i];
                    }
                    }

                    int driven=0;
                    double s=0.0;
                    int j=n;

                    if (!driven) 
                    {
                        s=b[index[n]]/delta[n];
                        if (flag) cDynPrintf("init s set to %5.9g delta=%5.9g\n",s,delta[n]);
                    }


                    // next check active constraints
                    if (s >= 0.0f) 
                    {
                        for (int i=n-1;i>=0;i--) 
                        {
                            if (BASE(type[index[i]][pass]) == CDYN_UNILATERAL && delta[i] > 0.0f) 
                            {
                                double sp = x[index[i]]/delta[i];
                                if (sp < s) 
                                {
                                    if (flag) cDynPrintf("change s set to %5.9g delta=%5.9g %d\n",s,delta[i],index[i]);
                                    s=sp;
                                    j=i;
                                }
                            }
                        }
                    } 
                    else 
                    { 
                        // need to handle BILATERAL CONSTRAINT
                        for (int i=n-1;i>=0;i--) 
                        {
                            if (BASE(type[index[i]][pass]) == CDYN_UNILATERAL && delta[i] < 0.0f) 
                            {
                                double sp = x[index[i]]/delta[i];
                                if (sp > s) 
                                {
                                    s=sp;
                                    j=i;
                                }
                            }
                        }
                    }

                    // next check if any constraints saturate
                    int saturation=0;
                    for (int i=n;i>=0;i--) 
                    {
                    //  if (type[index[i]][pass] == CDYN_FRICTION + CDYN_BILATERAL) {
                        if (BASE(type[index[i]][pass]) == CDYN_BILATERAL && max[index[i]] > 0.0) 
                        {
                            double deltafn=(i == n)?1:-delta[i];
                            if (fabs(deltafn) < 1e-6f) continue;
                            double sp=(((deltafn*s < 0.0f)?-max[index[i]]:max[index[i]])-x[index[i]])/deltafn;
                            
                            if (flag) 
                            {
#if 0
                                cDynPrintf("Delta\n");
                                for (j=0;j<n;j++) 
                                {
                                    cDynPrintf("%13d, ", index[i]);
                                }
                                cDynPrintf("\n");
                                for (j=0;j<n;j++) 
                                {
                                    cDynPrintf("%13.5g, \n", delta[i]);
                                }
#endif
                                cDynPrintf("%d s=%5.5f sp=%5.5f deltafn=%5.5g\n",index[i], s,sp, deltafn);
                            }
                            
                            if ((sp > 0.0 && sp < s) || (sp < 0.0 && sp > s)) 
                            {
/*
                                if (FLAG(type[index[i]][pass]) != CDYN_FRICTION)
                                    cDynPrintf("saturation(%5.17f)\n",max[index[i]]);
*/
                                saturation++;
                                s=sp;
                                j=i;
                            }
                        }
                    } // for saturation

#ifdef RTIMER1
                    rtimer.start();
#endif
                    if (flag) cDynPrintf("Block Variable %d\n",index[j]);

//---------------------------------------------------------------------------
#if 0
//---------------------------------------------------------------------------

                    double delta2[m];
                    {for (int i=0;i<m;i++) 
                    {
                        delta2[i]=delta[i];
                    }}

                    // update a different way
                    {for (int i=0;i<=n;i++) 
                    {
                        int x=index[i];
                        int xx=x;
                        while (contact[xx] == NULL) xx--;
                        cDynContactList* list;
                        for (int j=0;j<2;j++) 
                        {
                            if (j==0) list=contact[xx]->list_;
                            else if (contact[xx]->pair_ == NULL) continue;
                            else list=contact[xx]->pair_->list_;
                            for (int k=0;k<list->n_;k++) {
                                int y=list->list_[k]->n_;
                                if (list->list_[k]->prim_ != NULL)
                                    for (int l=0;l<list->list_[k]->numPoints_;l++) {
                                        if (rindex[y] > n)
                                            delta2[rindex[y]] += (i==n)?-lambda(x,y):lambda(x,y)*delta[i];
                                        y++;
                                    }
                            }
                        }
                    }}

//---------------------------------------------------------------------------
#endif
//---------------------------------------------------------------------------

                    // update delta for constraints not yet tested
                    {for (int k=n+1;k<m;k++) 
                    {
                        delta[k]= -lambda(k,n,index);
                    }}

                    {for (int i=0;i<n;i++) 
                    {
                        double d=delta[i];
                        for (int k=n+1;k<m;k++) 
                        {
                            delta[k] += lambda(k,i,index)*d;
                        }
                    }}
#if 0
                    // update delta for constraints not yet tested
                    {for (int k=n+1;k<m;k++) 
                    {
                        delta[k]= -lambda(k,n,index);
                        for (int i=0;i<n;i++) 
                        {
                            delta[k] += lambda(k,i,index)*delta[i];
                        }
                    }}
#endif
#if 0
                    // update delta for constraints not yet tested
                    {for (int k=n+1;k<m;k++) 
                    {
                        delta[k]= -lambda(k,n,index);
                        for (int i=0;i<n;i++) 
                        {
                            delta[k] += lambda(k,i,index)*delta[i];
                        }
                    }}
#endif

#ifdef RTIMER1
                    rtimer.end();
                    rtotal += rtimer.time();
                    rcount++;
#endif
                    if (j==n && s==0.0f && BASE(type[index[j]][pass]) == CDYN_BILATERAL) 
                    {
                        assert(b[index[j]]==0.0f && x[index[j]]==0.0f);
                        n++;
                        continue;
                    }

                    // NO LONGER applies assert(BASE(type[index[j]][pass])==CDYN_BILATERAL||s>=0.0f);

                    if (s==0.0f) 
                    {
                        b[index[j]]=0.0f; // has to be zero but sometimes we get numerical errors
                        assert(b[index[j]]==0.0f && x[index[j]]==0.0f);
                    } 
                    else 
                    { 
                        // s!=0

                        //-------------------------------------------------------------------
                        //
                        //  UPDATE x and b given s
                        //
                        //-------------------------------------------------------------------

                        // update x and b for active constraints
                        {for (int i=0;i<n;i++) 
                        {
                            x[index[i]] -= s*delta[i];
                            if (BASE(type[index[i]][pass]) != CDYN_BILATERAL)
                            {
                                //assert(x[index[i]]>=0.0f);
                                if (x[index[i]]<0.0f)
                                {
                                    //cDynPrintf("x=%13.9f\n",x[index[i]]);
                                    //comment out for -g
                                    //assert(x[index[i]]>-CDYN_ERROR_BOUND); 
                                    x[index[i]]=0.0f;
                                }
                            }
                            assert(b[index[i]]==0.0f);
                        }
                        }

                        // and the new constraint
                        x[index[n]] += s; // note + is correct because delta = 1 not -1 as with others
                        if (BASE(type[index[n]][pass]) != CDYN_BILATERAL && x[index[n]]<0.0f) 
                        {
                            cDynPrintf("x=%5.9f\n",x[index[n]]);
                        }
                        assert(BASE(type[index[n]][pass]) == CDYN_BILATERAL || x[index[n]]>=0.0f);

                        b[index[n]] -= s*delta[n];

                        // b[index[n]] is driven to be <= 0.0 but because of numerical
                        // errors sometimes it is not.
                        // An element in C was selected when the new element should
                        // have been (option force j=n and b[index[n]]=0)
                        //assert(b[index[n]]<CDYN_ERROR_BOUND);  CHECK XXXXX FIX
#if 1
                        // and the other constraints
                        {for (int i=n+1;i<m;i++) 
                        {
                            //assert(x[index[i]] == 0.0f);
                            b[index[i]] -= s*delta[i];
                        }
                        }
#endif
                    } // if (s != 0)

                    //-------------------------------------------------------------------
                    //
                    // UPDATE SET C and NC
                    //
                    //-------------------------------------------------------------------

                    if (driven && j==n) 
                    {
                        cDynPrintf("FORCE DRIVEN %d\n",index[n]);
                        type[index[n]][pass]=FLAG(type[index[n]][pass]) + CDYN_NOTACTIVE;
                        nv=n;
                        dindex=-1;
                    } 
                    else if (saturation && j == n) 
                    {
                        //type[index[n]][pass]=FLAG(type[index[n]][pass]) + CDYN_NOTACTIVE;
         
                        for (int p=pass; p < CDYN_ERROR; p++)
                        type[index[n]][p]=FLAG(type[index[n]][p]) + CDYN_NOTACTIVE;
         
                        nv=n;
                        dindex=-1;
                    } 
                    else if (j == n) 
                    {
                        // new constraint goes to C
                        // assert(fabs(b[index[n]])<CDYN_ERROR_BOUND); CHECK XXX FIX
                        b[index[n]]=0.0f;

                        assert(BASE(type[index[j]][pass])!=CDYN_NOTACTIVE);

                        dindex= -1; // finished driving contraint to zero find new constraint
                        n++;
                    } 
                    else 
                    { 
                        // j is in C goes to NC
                        assert(j < n);
 
                        assert(b[index[j]]==0.0f);

                        //assert(BASE(type[index[j]][pass])!=CDYN_BILATERAL);
                        if (saturation) 
                        {
                            type[index[j]][pass]=FLAG(type[index[j]][pass]) + CDYN_NOTACTIVE;
                        } 
                        else 
                        {
                            x[index[j]]=0.0f;
                        }

                        loopCount--;
                        int tmp=index[j];
                        dindex=n; // continue next loop with same constraint
                        n--;
#if 0
                        index[j]=index[n];
                        rindex[index[j]]=j;
#else
                        {for (int i=j;i<n;i++)
                        {
                            index[i] = index[i+1];
                            rindex[index[i]] = i;
                        }
                        }
#endif
                        index[n]=tmp;
                        rindex[index[n]]=n;
                        nv=j;

                        loopCount--;
                    } 
                }
                
                // BEGIN DEBUG
                if (flag) 
                {
                    cDynDebugResolve(m,n,index,type,b,x);
                } 
                // END DEBUG
            } //  if (mindex == -1) else 

#ifdef CHECKMATH
            // DEBUG CHECK math
            {
                static double maxdiff=0.0;
                for (int i=0;i<m;i++) 
                {
                    double xx=debug[i];
                    for (int j=0;j<m;j++) 
                    {
                        xx += lambda(i,j)*x[j];
                    }
                    double diff = fabs(b[i] - xx);
                    if (diff > maxdiff) 
                    {
                        maxdiff = diff;
                        cDynPrintf("MATH DIFF %5.9g\n",maxdiff);
                    }
                }
            }
            // END DEBUG
#endif
      
        } // while (!done)

        if (flag && loopCount > m)
        {
            cDynPrintf("  **** loopCount(%d)-m(%d)=%d ****  ",loopCount,m,loopCount-m);
        }

        // CHECK TO SEE IF ANY SLEEPING BASENODES WOULD HAVE BEEN AWAKEN BY ACTION
        bool awaken=false;
        if (pass == CDYN_IMPULSE) 
        {
            for (int i=0; i < m; ) {
                cDynContact* c=contact[i];
                currentNum=i;
                if (c->type_ == CDYN_CONTACT) 
                {
                    if (c->pair_ == NULL && c->prim_->nodeB_ != NULL && c->prim_->nodeB_->baseNode()->status() == CDYN_SLEEP) 
                    {
                        int ok=0;
                        double value=(pass == CDYN_IMPULSE)?c->prim_->nodeB_->baseNode()->contact()->wakeImpulse():c->prim_->nodeB_->baseNode()->contact()->wakeForce();
                        for (int j=0;j<c->numPoints_;j++) 
                        {
                            if (fabs(x[i+j]) > value) { ok++; break; }
                        }
                        if (ok) 
                        {
                            awaken=true;
                            c->prim_->nodeB_->baseNode()->status(CDYN_ACTIVE);
                            if (flag) {
                                cDynPrintf("BaseNode %s being set to ACTIVE status\n",
                                    c->prim_->nodeB_->baseNode()->data());
                            }
                        }
                    }
                    i += c->numPoints_;
                } 
                else 
                    i++;
            }
        }

        if (awaken) 
        {
            if (bigLambda) free(Lambda);
            return(false); // system has changed rerun with new ACTIVE base nodes
        }
      
        //-------------------------------------------------------------------
        //
        // FOUND SOLUTION, APPLY IMPULSE/FORCE/ERROR TO SYSTEM
        //
        //-------------------------------------------------------------------
        
        currentX=x;
        {for (int i=0; i < m; ) 
        {
            cDynContact* c=contact[i];
            currentNum=i;
            if (c->type_ == CDYN_CONTACT) 
            {
                switch (pass) 
                {
                    case CDYN_IMPULSE:
                    // Apply computed impulses back to system.
                    //		      c->list_->base_->fwdCallbacks(
                    //			&DeCallbackJtIn, &cDynCallbackJtSphereIn,
                    //			&DeCallbackVelocityOut, &DeCallbackVelocitySphereOut);
                    //		      if (c->pair_ != NULL) c->pair_->list_->base_->fwdCallbacks(
                    //			&DeCallbackJtIn, &cDynCallbackJtSphereIn,
                    //			&DeCallbackVelocityOut, &DeCallbackVelocitySphereOut);
                    break;
                    case CDYN_FORCE:
                    // Apply computed forces back to system.
                    //		      c->list_->base_->fwdCallbacks(
                    //			&DeCallbackJtIn, &cDynCallbackJtSphereIn,
                    //			&DeCallbackAccelerationOut, &DeCallbackAccelerationSphereOut);
                    //		      if (c->pair_ != NULL) c->pair_->list_->base_->fwdCallbacks(
                    //			&DeCallbackJtIn, &cDynCallbackJtSphereIn,
                    //			&DeCallbackAccelerationOut, &DeCallbackAccelerationSphereOut);
                    break;
                    case CDYN_ERROR:
                    // Apply computed errors back to system.
                    //		      c->list_->base_->fwdCallbacks(
                    //			&DeCallbackJtIn, &cDynCallbackJtSphereIn,
                    //			&DeCallbackPositionOut, &DeCallbackPositionSphereOut);
                    //		      if (c->pair_ != NULL) c->pair_->list_->base_->fwdCallbacks(
                    //			&DeCallbackJtIn, &cDynCallbackJtSphereIn,
                    //			&DeCallbackPositionOut, &DeCallbackPositionSphereOut);
                    break;
                }
          
                int ok=0;

                for (int j=0;j<c->numPoints_;j++) 
                    if (x[i+j] != 0.0f) { ok++; break; }

                if (ok) 
                {
                    //cDynPrintf("setting(%d)\n",n_);

                    if (c->pair_ != NULL && c->list_->base_ != c->pair_->list_->base_) 
                    {
                        // specify a dependency between the bases for rollback
                        // XXX FIX :should have a flag to avoid marking dependency for both impulse/force
                        c->list_->base_->dependent(c->pair_->list_->base_);
                    }

                    currentContactPointDirection= -1.0f;
                    if (pass == CDYN_FORCE)
                    cDynamicsFwdDynamicsIn(c->prim_->nodeA_, &cDynCallbackJtIn, &cDynCallbackJtSphereIn);
                    
                    else if (pass == CDYN_IMPULSE) 
                    {
                        cDynamicsFwdDynamicsConfig(c->list_->base_->root(),c->prim_->nodeA_, &cDynCallbackJtIn, &cDynCallbackJtSphereIn, &cDynCallbackVelocityOut, &cDynCallbackVelocitySphereOut, 1);
                    }
                    
                    else
                    {
                        cDynamicsFwdDynamicsConfig(c->list_->base_->root(),c->prim_->nodeA_, &cDynCallbackJtIn, &cDynCallbackJtSphereIn, &cDynCallbackPositionOut, &cDynCallbackPositionSphereOut, 0);
                    }
                    
                    if (c->pair_ != NULL) 
                    {
                        currentContactPointDirection=  1.0f;
                        if (pass == CDYN_FORCE)
                            cDynamicsFwdDynamicsIn(c->prim_->nodeB_, &cDynCallbackJtIn, &cDynCallbackJtSphereIn);
                        else if (pass == CDYN_IMPULSE)
                            cDynamicsFwdDynamicsConfig(c->pair_->list_->base_->root(),c->prim_->nodeB_, &cDynCallbackJtIn, &cDynCallbackJtSphereIn, &cDynCallbackVelocityOut, &cDynCallbackVelocitySphereOut, 1);
                        else
                            cDynamicsFwdDynamicsConfig(c->pair_->list_->base_->root(),c->prim_->nodeB_, &cDynCallbackJtIn, &cDynCallbackJtSphereIn, &cDynCallbackPositionOut, &cDynCallbackPositionSphereOut, 0);
                    }

                    if (pass == CDYN_ERROR) { // I believe these functions can be eliminated in the future FIX XXX
                    c->list_->error(true);
                    if (c->pair_ != NULL) c->pair_->list_->error(true);
                }
            }
            i += c->numPoints_;
        } 
        else 
        { 
            // joint limit
            switch (pass) 
            {
            case CDYN_IMPULSE:
            //	          c->list_->base_->fwdCallbacks(
            //		  	&DeCallbackJtIn, &cDynCallbackJtSphereIn,
            //		  	&DeCallbackVelocityOut, &DeCallbackVelocitySphereOut);
                if (x[i] != 0.0f)
                    cDynamicsFwdDynamicsConfig(c->list_->base_->root(), c->joint_->object(), &cDynCallbackJtIn, &cDynCallbackJtSphereIn, &cDynCallbackVelocityOut, &cDynCallbackVelocitySphereOut, 1);
                break;

            case CDYN_FORCE:
            //		  c->list_->base_->fwdCallbacks(
            //		  	&DeCallbackJtIn, &cDynCallbackJtSphereIn,
            //		  	&DeCallbackAccelerationOut, &DeCallbackAccelerationSphereOut);
                if (x[i] != 0.0f)
                    cDynamicsFwdDynamicsConfig(c->list_->base_->root(), c->joint_->object(), &cDynCallbackJtIn, &cDynCallbackJtSphereIn, &cDynCallbackAccelerationOut, &cDynCallbackAccelerationSphereOut, 1);
                break;

            case CDYN_ERROR:
            //		  c->list_->base_->fwdCallbacks(
            //		  	&DeCallbackJtIn, &cDynCallbackJtSphereIn,
            //		  	&DeCallbackPositionOut, &DeCallbackPositionSphereOut);
                if (x[i] != 0.0f) 
                {
                    cDynamicsFwdDynamicsConfig(c->list_->base_->root(), c->joint_->object(), &cDynCallbackJtIn, &cDynCallbackJtSphereIn, &cDynCallbackPositionOut, &cDynCallbackPositionSphereOut, 0);
                    c->list_->error(true); // FIX: I believe this flag can be eliminated
                }
                break;
            }
            i++;
        }
    }
    }

    //-----------------------------------------------------------------------
    // BEGIN DEBUG
    //-----------------------------------------------------------------------

    if (flag) 
    {
        switch (pass) 
        {
            case CDYN_IMPULSE: cDynPrintf("IMPULSE\n");break;
            case CDYN_FORCE:   cDynPrintf("FORCE\n");break;
            case CDYN_ERROR:   cDynPrintf("ERROR\n");break;
        }
        {for (int i=0;i<m;) 
        {
            if (contact[i]->type_ != CDYN_CONTACT) 
            {
                cDynPrintf("%s [%12.8f,%12.8f,%12.8f]\n",
                    contact[i]->joint_->data(),
                    contact[i]->joint_->q(),
                    contact[i]->joint_->v(),
                    contact[i]->joint_->a()
                    );
                i++;
            } 
            else 
            {
                i += contact[i]->numPoints_;
            }
        }
        }
    }

    //-----------------------------------------------------------------------
    // END DEBUG
    //-----------------------------------------------------------------------
    
    // prepare for next loop
    switch (pass) 
    {
        case CDYN_IMPULSE: // prepare for force pass
#ifdef C_DYN_DEBUG_OPS
        {
            // compute the change in kinetic energy
            double vm[n];
            {for (int i=0;i<n;i++) 
            {
                vm[i]=vminus[index[i]];
            }}

            // forward subsitution
            {for (int j=0;j<n;j++) 
            {
                double ll=vm[j];
                for (int i=0;i<j;i++)
                    ll -= L(j,i,index)*vm[i];
                vm[j]= ll/diag[index[j]];
            }
            }
            // backward subsitution
            {for (int j=n-1;j>=0;j--) 
            {
                double ll=vm[j];
                for (int i=j+1;i<n;i++)
                    ll -= L(i,j,index)*vm[i];
                vm[j] = ll/diag[index[j]];
            }
            }

            double DeltaKE=0.0f;
            {for (int i=0;i<n;i++) 
            {
                DeltaKE += vm[i]*vminus[index[i]];
            }}

            DeltaEnergy += 0.5f*DeltaKE;
        }
#endif
        // find updated velocity
        {for (int i=0;i<m;i++) 
        {
            impulse[i]=x[i];
            v[i] += b[i];
            dv[i] = b[i];
        }}

        // zero the previous constraint torques: XXX FIX need to optimized
        // XXX: maybe redundant
        {for (int i=0;i<m;) 
        {
            cDynContact* c=contact[i];
            if (c->type_ == CDYN_CONTACT) 
            {
                if (!c->prim_->nodeA_->isFixed()) cDynamicsTorqueExternalZeroInwardPath(c->prim_->nodeA_);
                if (!c->prim_->nodeB_->isFixed()) cDynamicsTorqueExternalZeroInwardPath(c->prim_->nodeB_);
            // FIX set force to zero
            i += c->numPoints_;
            } 
            else 
            {
                c->joint_->torqueExternalZero();
                c->tau_=0.0f;
                i++;
            }
        }
        }
          
        // compute contact space acceleration
        acceleration(m,contact,b);
        {for (int i=0;i<m;i++) 
        {	
            x[i]=0.0f;
        }}
          
        // TEST HOO for dynamic friction
/*
        {for (int i=0;i<m;i++) 
        {
            if (type[i][CDYN_FORCE] == CDYN_FRICTION + CDYN_NOTACTIVE) 
            {
            b[i] = 0.1*v[i];
            }
        }}
*/
        break;

        case CDYN_FORCE: 
            // prepare for ERROR pass
            {for (int i=0;i<m;i++) 
            {
                // save acceleration and contact forces
                a[i]=b[i];	
                force[i]=x[i];
                b[i]=err[i];	// set b to error vector
                x[i]=0.0f;
            }}
            break;

        case CDYN_ERROR: break; // do nothing
    }
    } // multi-pass


    //---------------------------------------------------------------------------
    //
    // FIND MAXIMUM TIME STEP FOR NEXT INTEGRATION STEP
    //
    //---------------------------------------------------------------------------

    {for (int i=0;i<m;) 
    {
        cDynContact* c=contact[i];
        if (c->type_ == CDYN_CONTACT) 
        {
            cDynContactPoint* p=c->head_;
            while (p != NULL) 
            {
                if (BASE(type[i][CDYN_IMPULSE]) == CDYN_UNILATERAL && a[i] < -1e-8f) 
                {
                    // use quadatic rule to find maximum time before
                    // constraint would pass -0.5*maxerr[i]
                    if (err[i] < -maxerr[i]) 
                        err[i]= -maxerr[i];
                    
                    //double tmax= (-v[i]-cDynSqrt(v[i]*v[i]-2.0f*a[i]*(err[i]+maxerr[i])))/a[i];
                    //double tmax=(v[i]-cDynSqrt(v[i]*v[i]-a[i]*(b[i]+maxerr[i])))/a[i];
                    //double tmax= -v[i]/a[i];
                    // try find time when constraint will be satified exactly
                    
                    double sq=v[i]*v[i]-2.0f*a[i]*b[i];
                    double tmax=0.0f;
                    if (sq >= 0.0f) 
                    { 
                        // found time when constraint will be satisfied
                        tmax= (-v[i]-cDynSqrt(sq))/a[i];
                    } 
                    else if (v[i] >= 0.0f) 
                    { 
                        // if it never reaches position determine when velocity will be equal to zero
                        tmax= -v[i]/a[i];
                    } 
                    else 
                    { 
                        //if velocity is always negative, determine when maxerr is achieved
                        tmax= (-v[i]-cDynSqrt(v[i]*v[i]-2.0f*a[i]*(b[i]+maxerr[i])))/a[i];
                    }
                    
                    if (tmax < c->list_->inc_) 
                    {
                        c->list_->inc_=tmax;
                        //if (flag) cDynPrintf("STEP time=%-13.9f %s<--->%s (0x%08xl)\n", c->list_->time_+c->list_->inc_,
                        //	c->prim_->nodeA_->data(),c->prim_->nodeB_->data(),(unsigned long)c->prim_); 
                        if (flag) cDynPrintf("STEP time=%-13.9f %s<--->%s)\n", c->list_->time_+c->list_->inc_,
                        cDynObjectName(c->prim_->nodeA_),cDynObjectName(c->prim_->nodeB_)); 
                    }
                    if (c->pair_ != NULL && tmax < c->pair_->list_->inc_) 
                    {
                        c->pair_->list_->inc_=tmax;
                        //if (flag) cDynPrintf("STEP time=%-13.9f %s<--->%s (0x%08xl)\n", c->list_->time_+c->list_->inc_,
                        //	c->prim_->nodeA_->data(),c->prim_->nodeB_->data(),(unsigned long)c->prim_); 
                        if (flag) cDynPrintf("STEP time=%-13.9f %s<--->%s\n", c->list_->time_+c->list_->inc_,
                        cDynObjectName(c->prim_->nodeA_),cDynObjectName(c->prim_->nodeB_)); 
                    }
                }
                p=p->next_;i++;
            }
        } 
        else 
        {
            if (BASE(type[i][CDYN_IMPULSE]) == CDYN_UNILATERAL && a[i] < -1e-8) 
            {
                if (err[i] < -maxerr[i]) err[i]= -maxerr[i];

                //double tmax= (-v[i]-cDynSqrt(v[i]*v[i]-2.0f*a[i]*(err[i]+maxerr[i])))/a[i];
                //double tmax= -v[i]/a[i];
                // try find time when constraint will be satified exactly
                
                double sq=v[i]*v[i]-2.0f*a[i]*b[i];
                double tmax=0.0f;
                
                if (sq >= 0.0f) 
                { 
                    // found time when constraint will be satisfied
                    tmax= (-v[i]-cDynSqrt(sq))/a[i];
                } 
                else if (v[i] >= 0.0f) 
                { 
                    // if it never reaches position determine when velocity will be equal to zero
                    tmax= -v[i]/a[i];
                } 
                else 
                { 
                    //if velocity is always negative, determine when maxerr is achieved
                    tmax= (-v[i]-cDynSqrt(v[i]*v[i]-2.0f*a[i]*(b[i]+maxerr[i])))/a[i];
                }

                if (tmax < c->list_->inc_) c->list_->inc_=tmax;
                if (flag) cDynPrintf("STEP_J time=%-13.9f %s\n", c->list_->time_+c->list_->inc_,
                                c->joint_->data()); 
          }
          i++;
        }
    }}

    if (false) 
    {
        cDynPrintf(" impulse=%d force=%d error=%d cN=%d cZero=%d cNV=%d\n",cImpulse,cForce,cError,cN,cZero,cNV);
        cDynPrintf("The END!\n");
        //base_->root()->displaySmall();
        // compute contact space acceleration to see if worked it out correctly
        static double DEBerr=0.0f;
        setupReset();

        acceleration(m,contact,b);
        // print stuff
        cDynPrintf("DEBX  ");
        {for (int i=0;i<m;i++) 
        {
                char c=(i<n)?'C':'_';
                char t='E';
                switch (BASE(type[index[i]][pass])) 
                {
                    case CDYN_UNILATERAL: t= ' '; break;
                    case CDYN_BILATERAL:  t= '*'; break;
                    case CDYN_REDUNDANT:  t= 'R'; break;
                    case CDYN_NOTACTIVE:  t= '#'; break;
                }
                char f=' ';
                switch (FLAG(type[index[i]][pass])) 
                {
                    case CDYN_FRICTION_NORMAL: f= '^'; break;
                    case CDYN_SATURATED_NORMAL: f= 'S'; break;
                    case CDYN_FRICTION: f= 'F'; break;
                }
                cDynPrintf("%10d%c%c%c ",index[i], c, f, t);
        }
        }

        cDynPrintf("\na     ");
        {for (int i=0;i<m;i++) 
        {
            cDynPrintf("%13.9f ", b[index[i]]);
        }
        }
        cDynPrintf("\nA     ");
        {for (int i=0;i<m;i++) 
        {
            cDynPrintf("%13.9f ", a[index[i]]);
        }
        }
        cDynPrintf("\nx     ");
        {for (int i=0;i<m;i++) 
        {
            cDynPrintf("%13.9f ", force[index[i]]);
        }
        }
        cDynPrintf("\n");
        {for (int i=0;i<m;i++) 
        {
            if (fabs(b[index[i]] - a[index[i]]) > DEBerr) 
            {
                DEBerr=fabs(b[index[i]] - a[index[i]]);
                cDynPrintf("MAX DIFF %5.9f\n",DEBerr);
            }
        }
        }
    }
    
    // update f-curves
    update();

    //---------------------------------------------------------------------------
    //
    // REMOVE CONSTRAINTS THAT ARE NO LONGER ACTIVE
    //
    //---------------------------------------------------------------------------

    //if (false) // XXX: fix
    {for (int i=0;i<m;) 
    {
        cDynContact* c=contact[i];
        if (c->type_ == CDYN_CONTACT) 
        {
            cDynContactPoint* p=c->head_;
            while (p != NULL) 
            {
                cDynContactPoint* t=p;
                p=p->next_;
                // XXX: maybe residue left over
                if (BASE(type[i][CDYN_IMPULSE]) == CDYN_REDUNDANT || ( fabs(force[i]) < CDYN_ERROR_BOUND && fabs(impulse[i]) < CDYN_ERROR_BOUND)) 
            { 
                // delete a contact point
                delete t;
            } 
                else 
            { 
                // update the list
                t->impulse_ = impulse[i];
                t->force_ = force[i];
                t->velocity_ = v[i];
                t->deltavelocity_ = dv[i];
                t->acceleration_ = a[i];
            }
                i++;
            }
            if (c->numPoints_ == 0) 
            { 
                // if no points left get rid of it
                c->list_->remove(c->index_);
            }
        } 
        else 
        { 
            // joint limit
            // if (v[i] >= CDYN_ERROR_BOUND || fabs(force[i]) < CDYN_ERROR_BOUND) 
            if (BASE(type[i][CDYN_IMPULSE]) == CDYN_REDUNDANT || fabs(force[i]) < CDYN_ERROR_BOUND) 
            {
                c->list_->remove(c->index_);
            } 
            else 
            {
                c->tau_ = force[i];
                //c->joint_->a(0.0f);
                //c->joint_->v(0.0f);
                //c->joint_->status((c->type_ == CDYN_JOINT_LOWER)?CDYN_AT_LOWER:CDYN_AT_UPPER);
            }
            i++;
        }
    }
    }

#if 0
    {for(int y=0;y<m;y++)
    {
        if (BASE(type[y][CDYN_FORCE]) != CDYN_REDUNDANT && BASE(type[y][CDYN_FORCE]) != CDYN_NOTACTIVE) {
            if (!((a[y]>= -1e-5 || BASE(type[y][CDYN_FORCE])==CDYN_NOTACTIVE) && (force[y]>=0.0f || BASE(type[y][CDYN_FORCE])==CDYN_BILATERAL || BASE(type[y][CDYN_FORCE])==CDYN_NOTACTIVE)))
        {
            assert(a[y]*force[y]==0.0f);
            if (BASE(type[y][CDYN_FORCE])==CDYN_BILATERAL)
            cDynPrintf("CDYN_BILATERAL\n");
            else
            cDynPrintf("NOT CDYN_BILATERAL(A) -- %d\n",BASE(type[index[y]][CDYN_FORCE]));
            cDynPrintf("%d, a=%13.9f, f=%13.9f\n",y,a[y],force[y]);
            exit(1);
        }
        if (BASE(type[y][CDYN_IMPULSE]) != CDYN_REDUNDANT && BASE(type[y][CDYN_IMPULSE]) != CDYN_NOTACTIVE) {
            assert(v[y] >= -1e-5);
        }
        }
    }
    }
#endif

    if (bigLambda) free(Lambda);
    //cDynPrintf("setting finish(%d)\n",n_);

    return(true);
}

//---------------------------------------------------------------------------

void cDynContactList::display(const int m, const int n, double* l, double *diag, double *b, int* index)
{
    {for (int i=0;i<n;i++) 
    {
        for (int j=0;j<n;j++) 
        {
            cDynPrintf("%6.3f ",lambda(i,j,index));
        }
        cDynPrintf("      %6.3f\n",b[index[i]]);
    }}
    cDynPrintf("\n");
    if (diag != NULL) 
    {
        {for (int i=0;i<n;i++) 
        {
                cDynPrintf("%6.3f ",diag[index[i]]);
        }}
    }
    cDynPrintf("\n\nindex ");
    {for (int i=0;i<n;i++) 
    {
            cDynPrintf("%d ",index[i]);
    }}
    cDynPrintf("\n");
}

//---------------------------------------------------------------------------

void cDynContactList::display(const int m, double* l, double *b)
{
    for (int k=0;k<m;k += NCOL) 
    {
        if (m > NCOL) 
        {
            cDynPrintf("Cols: %d-%d\n",k,(k+NCOL < m)?k+NCOL-1:m);
        }
        
        for (int i=0;i<m;i++) 
        {
            int j=k;
            for (j=k;j<m && j<k+NCOL;j++) 
            {
                cDynPrintf("%12.9f ",lambda(i,j));
            }
            if (j == m)
                cDynPrintf("      %12.9f\n",b[i]);
            else
                cDynPrintf("\n");
        }
    }
}

//---------------------------------------------------------------------------
#define SNCOL	256
//---------------------------------------------------------------------------

void cDynContactList::displaySmall(const int m, double* l, double *b)
{
    for (int k=0;k<m;k += SNCOL) 
    {
        if (m > SNCOL)
        {
            cDynPrintf("Cols: %d-%d\n",k,(k+SNCOL < m)?k+SNCOL-1:m);
        }
        for (int i=0;i<m;i++) 
        {
            int j=k;
            for (j=k;j<m && j<k+SNCOL;j++) 
            {
                cDynPrintf("%c",(lambda(i,j)==0.0)?' ':'X');
            }
            cDynPrintf("\n");
        }
    }
}

//---------------------------------------------------------------------------


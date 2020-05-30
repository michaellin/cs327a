//===========================================================================
/*
    This file is part of the dynamics library.
    Copyright (C) 2014, Artificial Intelligence Laboratory,
    Stanford University. All rights reserved.

    \version   $MAJOR.$MINOR.$RELEASE $Rev: 1242 $
*/
//===========================================================================

//---------------------------------------------------------------------------
#include <stdlib.h>
#include <string.h>
#include "utility/CDynLogger.h"
#include "matrix/CDynTransform.h"
#include "object/CDynObject.h"
#include "node/CDynCollisionCheckRecord.h"
#include "node/CDynPrimPair.h"
#include "node/CDynFrictionRecord.h"
#include "node/CDynWorld.h"
//---------------------------------------------------------------------------
long cDynPrimPair::CALLNUM=1;
cDynPrimPair* cDynPrimPair::free_=NULL;
//---------------------------------------------------------------------------
// number used for hash table
#define DEPRIM_CACHESIZE 2749

// numbers should be relatively prime to each other. In this case all numbers are primes
#define A_MULTIPLIER 151
#define B_MULTIPLIER 1109
#define C_MULTIPLIER 751
//---------------------------------------------------------------------------
static const cDynTime tlimit=1e-8;
//---------------------------------------------------------------------------

class cDynProjectionVector 
{
    //-----------------------------------------------------------------------
    // CONSTRUCTOR & DESTRUCTOR:
    //-----------------------------------------------------------------------
public:

    //! Constructor of cDynProjectionVector. 
    cDynProjectionVector() {}

    //! Constructor of cDynProjectionVector. 
    cDynProjectionVector(const double x, const double y, const double err)
    { 
        data_[0]=x; 
        data_[1]=y; 
        data_[2]=err; 
    }


    //-----------------------------------------------------------------------
    // OPERATORS:
    //-----------------------------------------------------------------------
public:

    double &operator[](int i) { return data_[i]; }
    const double &operator[](int i) const { return data_[i]; }


    //-----------------------------------------------------------------------
    // PUBLIC MEMBERS:
    //-----------------------------------------------------------------------
public:

    void values(const double x, const double y, const double err) { data_[0]=x; data_[1]=y; data_[2]=err; }
    bool rightOf(const cDynProjectionVector& a, cDynProjectionVector& b) 
    {
        double xA=a[0]-data_[0]; double yA=a[1]-data_[1];
        double xB=b[0]-data_[0]; double yB=b[1]-data_[1];
        return((xA*yB - xB*yA) > 0.0f);
    }
    
    //---------------------------------------------------------------------------

    bool above(const cDynProjectionVector& a) 
    {
        return(a[0] == data_[0] && a[1] <= data_[1]);
    }
    
    //---------------------------------------------------------------------------

    bool below(const cDynProjectionVector& a) 
    {
        return(a[0] == data_[0] && a[1] >= data_[1]);
    }

    //---------------------------------------------------------------------------

    bool intersect(const cDynProjectionVector &pa, 
        const cDynProjectionVector &pb, 
        const cDynProjectionVector &qa, 
        const cDynProjectionVector& qb, 
        const int dir)
    {
        double dp[3]={pb[0]-pa[0],pb[1]-pa[1]};
        double dq[3]={qa[0]-qb[0],qa[1]-qb[1]};
        double dpq[2]={qa[0]-pa[0],qa[1]-pa[1]};
        double s= (dp[0]*dq[1] - dp[1]*dq[0]);
        if (fabs(s) < 1e-8) return(false);
        double u=(dq[1]*dpq[0] - dq[0]*dpq[1])/s;
        if (u >= 0.0f && u <= 1.0f) 
        {
            double v=(-dp[1]*dpq[0] + dp[0]*dpq[1])/s;
            if (v >= 0.0f && v <= 1.0f) 
            {
                dp[2]=pb[2]-pa[2];
                dq[2]=qb[2]-qa[2];
                data_[0]=pa[0]+u*dp[0];
                data_[1]=pa[1]+u*dp[1];
                double err=(qa[2]+v*dq[2]) - (pa[2]+u*dp[2]);
                data_[2]=dir?-err:err;
                return(true);
            }
        }
        return(false);
    }
    //---------------------------------------------------------------------------

    bool between(const cDynProjectionVector& n, 
        const cDynProjectionVector &pa,
        const cDynProjectionVector &pb, 
        const cDynProjectionVector &qa, 
        const cDynProjectionVector& qb, 
        const int dir) 
    {
        double dpn[2]={n[0]-pa[0],n[1]-pa[1]};
        double dp[3]={pb[0]-pa[0],pb[1]-pa[1],0.0f};
        double dqn[2]={n[0]-qa[0],n[1]-qa[1]};
        double dq[3]={qb[0]-qa[0],qb[1]-qa[1],0.0f};
        double up,uq;
        if (dp[0] < 1e-8f) up=(dp[1] > 0.0f)?1.0f:0.0f;
        else up=dpn[0]/dp[0];
        if (dq[0] < 1e-8f) uq=(dq[1] > 0.0f)?0.0f:1.0f;
        else uq=dqn[0]/dq[0];

        double yp=pa[1]+up*dp[1];
        double yq=qa[1]+uq*dq[1];

        if (yp <= n[1] || yq >= n[1]) return(false);

        data_[0]=n[0];
        data_[1]=n[1];

        dp[2]=pb[2]-pa[2];
        dq[2]=qb[2]-qa[2];

        double dy=yp-yq;
        double ep=pa[2]+up*dp[2];
        double eq=qa[2]+uq*dq[2];
        double v;
        if (dy < 1e-8f) v=0.5f;
        else v=(n[1] - yq)/dy;

        double err=eq + v*(ep - eq) - n[2];
        data_[2]=dir?-err:err;
        return(true);
    }


    //-----------------------------------------------------------------------
    // PRIVATE MEMBERS:
    //-----------------------------------------------------------------------
private:

    double data_[3];
};

//---------------------------------------------------------------------------

static cDynPrimPair cache[DEPRIM_CACHESIZE];

//---------------------------------------------------------------------------

static int cDynXCompare(const void* a, const void* b)
{
    double pA= (**(cDynProjectionVector **)a)[0];
    double pB= (**(cDynProjectionVector **)b)[0];
    //cDynPrintf("test %5.5f -- %5.5f\n",pA,pB);
    
    if (pA > pB) return(1);
    if (pA < pB) return(-1);
    return(0);
}

//---------------------------------------------------------------------------

static void cDynHull(cDynProjectionVector** list, 
    int n, 
    cDynProjectionVector** lower, 
    int& nl, 
    cDynProjectionVector** upper, 
    int& nu)
{
    int i; //index in sweep

    // add first two points (if they exist) to hull
    for (i=0;i<2 && i<n;i++) 
        lower[i]=upper[i]=list[i];

    nl=nu=i-1;
    for ( ;i<n;i++) 
    { 
        // insert remaining points updating accordingly
        
        // update upper bound
        while (nu > 0 && (list[i]->above(*upper[nu]) || upper[nu-1]->rightOf(*upper[nu],*list[i])))
            nu--;
        upper[++nu]=list[i];

        // update lower bound
        while (nl > 0 && (list[i]->below(*lower[nl]) || lower[nl-1]->rightOf(*list[i],*lower[nl])))
            nl--;
        lower[++nl]=list[i];
    }

    nl++;
    nu++;
}

//---------------------------------------------------------------------------

cDynPrimPair* cDynPrimPair::pairCache() const 
{
    unsigned long h=((
        A_MULTIPLIER*(nodeA_->uid()) ^ B_MULTIPLIER*((unsigned long)primA_->id) ^
        C_MULTIPLIER*(nodeB_->uid()) ^ ((unsigned long)primB_->id)
    )%DEPRIM_CACHESIZE);
    return(&(cache[h]));
}
    
//---------------------------------------------------------------------------

int cDynPrimPair::checkPrim(const cDynTransform& Tab) 
{
    static int iter=0;
    static cDynPrim* BinA=NULL;
    static int max=0;

    if (max < primB_->num) 
    {
        if (BinA != NULL) 
        {
            BinA->num=max;
            delete BinA;
            max=primB_->num;
        } 
        else 
        {
            max=256;
            if (max < primB_->num) max=primB_->num;
        }
        BinA=new cDynPrim(max);
    }

    // fetch witness information if available
    
    iter++;
    cDynPrimPair* entry=pairCache();
    if (entry->nodeA_ == nodeA_ && entry->nodeB_ == nodeB_ && entry->primA_ == primA_ && entry->primB_ == primB_) 
    {
        if (entry->callnum_ == CALLNUM) return(2); // already checked pair
        //CDYN_DISTW_NOV(&entry->witness_)=0;
    } 
    else 
    { 
        //reset witness
        entry->nodeA_=nodeA_;
        entry->primA_=primA_;
        entry->nodeB_=nodeB_;
        entry->primB_=primB_;
        CDYN_DISTW_NOV(&entry->witness_)=0;
    }
    entry->callnum_=CALLNUM;

    // rotate points in B to frame of A. Store results in BinA
    for (int i=0;i<primB_->num;i++) 
    {
        BinA->v[i].multiply(Tab, primB_->v[i]);
    }
    BinA->id=primB_->id;
    BinA->num=primB_->num;
    double r = primA_->r + primB_->r;
    double err = primA_->err + primB_->err;

//	if (cDynDebugLevel == 0)
//	{
//		cDynPrintf("A(%d,%d):%s\n",entry->primA_->id,entry->primA_->num,entry->primA_->id,entry->nodeA_->name()->string());
//		cDynPrintf("B(%d,%d):%s\n",entry->primB_->id,entry->primB_->num,entry->primB_->id,entry->nodeB_->name()->string());
//	}

//	dis_ = cDynDistDistance(primA_->v,primA_->num,BinA->v,BinA->num,&entry->witness_,r,true);
    dis_ = cDynDistDistance2(primA_->v,primA_->num,BinA->v,BinA->num,&entry->witness_,r,true);
    dis_ -= r;
    witness_= entry->witness_;
    
    if (dis_ > 0.0f) return(1);
    if (dis_ <= -err) return(-1);
    return(0);
}

//---------------------------------------------------------------------------

cDynTime cDynPrimPair::checkPrimTime(const cDynTime tl, 
    const cDynTime tu, 
    cDynTransform& Tab, 
    cDynPrimPairError& error)
{
    cDynTime tmin=tl;
    cDynTime tmax=tu;

    cDynFrame F;
    //cDynPrimPairCallnum();
    int type=checkPrim(Tab);
    if (type != 1) 
    { 
        // type ==1 means no penetration
        // otherwise use bisection to find contact time
        while (type != 0) 
        {
            cDynTime time=tmin + 0.5*(tmax - tmin);
            F.inversedMultiply(nodeA_->globalFrame(time), nodeB_->globalFrame(time));
            Tab.set(F);
            cDynPrimPairCallnum();
            type=checkPrim(Tab);

            if (type == 0) return(time);

            if (tmax - tmin < tlimit) 
            {
                if (tmin == tl) 
                {
                    // do a sanity check at time tl to make sure no penetration existed at that time
                    F.inversedMultiply(nodeA_->globalFrame(tl), nodeB_->globalFrame(tl));
                    Tab.set(F);
                    cDynPrimPairCallnum();
                    int typelow=checkPrim(Tab);

                    // if test fails then there is something definitly wrong
                    if (typelow == -1) { error = CDYN_UNKNOWN_FAILURE; return(tu); }

                    // if original range was to small report iteration regresion failure
                    if (tu - tl < tlimit) { error = CDYN_REGRESSION_FAILURE; return(tu); }
#if 0
                    // otherwise try droping the integration rate to see if problem can be resolved
                    cDynTime inc=0.5*(tu - tl);
                    cDynTime tn=tl+inc;
                    nodeA_->baseNode()->backup(tl);
                    nodeA_->baseNode()->integrate(tl,inc);
                    if (nodeA_->baseNode() != nodeB_->baseNode()) 
                    {
                        nodeB_->baseNode()->backup(tl);
                        nodeB_->baseNode()->integrate(tl,inc);
                    }
                    
                    cDynPrintf("ITERATION REGRESSION CONTACT %5.5f\n",inc);
                    F.inversedMultiply(nodeA_->globalFrame(tn), nodeB_->globalFrame(tn));
                    Tab.set(F);
                    tn=checkPrimTime(tl,tn,Tab,error);
                    return(tn);
#endif
                    cDynPrintf("ITERATION REGRESSION CONTACT\n");
                    error = CDYN_REGRESSION_FAILURE;
                    return(tu);
                } 
                else 
                { 
                    // The tolerances are to tight to find a contact time (for now mark as invalid)
                    error = CDYN_TOLERANCE_FAILURE; return(tu);
                }
            }

            if (type ==  1)  tmin=time;
            else tmax=time;
        }
    }
    return(tmax);
}

//---------------------------------------------------------------------------

int cDynPrimPair::findPatch(const cDynVector3& x, 
    const cDynVector3& n, 
    const double min, 
    const double max, 
    const double maxerr, 
    const double rA, 
    const double rB, 
    const cDynFrame& gA, 
    const cDynFrame& gB, 
    cDynContact* contact, 
    cDynFrictionRecord* fr)
{
    double d=n.dot(x);
    double dA=d + min;
    double dB=d + max;
    static int flag=0;

    int nA=0;
    int nB=0;
    cDynProjectionVector* pA=(cDynProjectionVector *)CDYN_ALLOCA(primA_->num*sizeof(cDynProjectionVector)); 
    cDynProjectionVector* pB=(cDynProjectionVector *)CDYN_ALLOCA(primB_->num*sizeof(cDynProjectionVector)); 
    cDynProjectionVector** ptrA=(cDynProjectionVector **)CDYN_ALLOCA(primA_->num*sizeof(cDynProjectionVector *)); 
    cDynProjectionVector** ptrB=(cDynProjectionVector **)CDYN_ALLOCA(primB_->num*sizeof(cDynProjectionVector *)); 

    // determine which axis to project points along
    // largest component goes to index[2]
    int index[3];
    if (fabs(n[0]) > fabs(n[1]))
        if (fabs(n[0]) > fabs(n[2])) { index[0]=1; index[1]=2;index[2]=0; }
        else {index[0]=0; index[1]=1; index[2]=2; }
    else 
        if (fabs(n[1]) > fabs(n[2])) { index[0]=0; index[1]=2; index[2]=1; }
        else {index[0]=0; index[1]=1; index[2]=2; }

    cDynVector3 pt;
    double err;
    double offset1=n[index[0]]*rA;
    double offset2=n[index[1]]*rA;
    
    for (int i=0;i<primA_->num;i++) 
    {
        pt.multiply(gA, primA_->v[i]);
        err=n.dot(pt) + rA;
        if (err >= dA) 
        {
            pA[nA].values(pt[index[0]] + offset1, pt[index[1]] + offset2, err - d);
            if (flag) cDynPrintf("%5.9f %5.9f %5.9f\n",pA[nA][0],pA[nA][1],pA[nA][2]);
            ptrA[nA]= &pA[nA]; nA++;
        }
    }

    offset1= -n[index[0]]*rB;
    offset2= -n[index[1]]*rB;
    
    for (int i=0;i<primB_->num;i++) 
    {
        pt.multiply(gB, primB_->v[i]);
        err=n.dot(pt) - rB;
        if (err <= dB) {
            pB[nB].values(pt[index[0]] + offset1, pt[index[1]] + offset2, err - d );
            if (flag) cDynPrintf("%5.9f %5.9f %5.9f\n",pB[nB][0],pB[nB][1],pB[nB][2]);
            ptrB[nB]= &pB[nB]; nB++;
        }
    }

    if (nA <= 1 || nB <= 1) return(0);

    // sort kists along x axis
    qsort(ptrA, nA, sizeof(cDynProjectionVector *), cDynXCompare);
    qsort(ptrB, nB, sizeof(cDynProjectionVector *), cDynXCompare);

    int boundN[2][2]; // number of vertices in envelope

    // boundN[OBJ[ENV]
    //   OBJ=0   primitive A
    //   OBJ=1   primitive B
    //   
    //   ENV=0   lower envelope
    //   ENV=1   upper envelope 

    // construct convex hulls
    cDynProjectionVector **boundA[2]=
    {
        (cDynProjectionVector **)CDYN_ALLOCA(nA*sizeof(cDynProjectionVector *)),
        (cDynProjectionVector **)CDYN_ALLOCA(nA*sizeof(cDynProjectionVector *))
    };

    cDynHull(ptrA, nA, boundA[0], boundN[0][0], boundA[1], boundN[0][1]);

    cDynProjectionVector **boundB[2]=
    {
        (cDynProjectionVector **)CDYN_ALLOCA(nB*sizeof(cDynProjectionVector *)),
        (cDynProjectionVector **)CDYN_ALLOCA(nB*sizeof(cDynProjectionVector *))
    };

    cDynHull(ptrB, nB, boundB[0], boundN[1][0], boundB[1], boundN[1][1]);

    cDynProjectionVector **bound[2][2]={{boundA[0],boundA[1]}, {boundB[0], boundB[1]}};

    // bound[OBJ][ENV][N] is a array of points
    //   OBJ=0   primitive A
    //   OBJ=1   primitive B
    //
    //   ENV=0   lower envelope
    //   ENV=1   upper envelope
    //   N= index number 0 >= N >= boundN[OBJ][ENV]

    int pos[2][2]={{0,0},{0,0}};
    
    // find intersection of the two hulls
    int pLeft=((*bound[0][0][0])[0] > (*bound[1][0][0])[0])?1:0;
    int pRight= pLeft?0:1;

    double xc=(*bound[pRight][0][0])[0];
    for (int i=0;i<2;i++) 
    {
        int j=1;
        while (j < boundN[pLeft][i] && (*bound[pLeft][i][j])[0] <= xc)
            j++;
        if (j == boundN[pLeft][i]) return(0); // no intersection
        pos[pLeft][i]=j;
        // assert(pos[pLeft][i] > 0);
    }

    cDynProjectionVector* clist=(cDynProjectionVector *)CDYN_ALLOCA(2*(nA+nB)*sizeof(cDynProjectionVector));

    int cnum=0;
    int c=pRight;
    int nc=pLeft;
    int side[2]={0,0}; // side[x]=y indicates for hull (x=0 A), (x=1 B)
               // which envelope (y=0 LOWER,y=1 UPPER) is to be
               // consider next (has smaller x value at pos[x][y]
    bool notline[2];
    notline[0]=(nA != 2);
    notline[1]=(nB != 2);

    pos[c][1]++; // do this to avoid counting first point twice

    // determine which of the envelopes is on the left
    side[nc]=((*bound[nc][0][pos[nc][0]])[0] < (*bound[nc][1][pos[nc][1]])[0])?0:1;
    for (;;) 
    {
        int e=side[c];
        // determine if new edge intersects edges of other hull

        if (pos[c][e] > 0) 
        {
            if (notline[nc] && clist[cnum].intersect(
                *bound[c][e][pos[c][e]-1], *bound[c][e][pos[c][e]],
                *bound[nc][0][pos[nc][0]], *bound[nc][0][pos[nc][0]-1],
                c
            )) cnum++;
            if (clist[cnum].intersect(
                *bound[c][e][pos[c][e]-1], *bound[c][e][pos[c][e]],
                *bound[nc][1][pos[nc][1]], *bound[nc][1][pos[nc][1]-1],
                c
            )) cnum++;
        }

        // determine if new point is between edges of other hull
        if (clist[cnum].between(*bound[c][e][pos[c][e]],
                        *bound[nc][1][pos[nc][1]-1],*bound[nc][1][pos[nc][1]],
                        *bound[nc][0][pos[nc][0]-1],*bound[nc][0][pos[nc][0]],
                        c
                )) cnum++;

        pos[c][e]++;

        if (pos[c][e] == boundN[c][e]) 
        {
            if (!notline[c]) break;
            e= e?0:1;
            while (pos[c][e] < boundN[c][e]) 
            {
                if (notline[nc] && clist[cnum].intersect(
                    *bound[c][e][pos[c][e]-1], *bound[c][e][pos[c][e]],
                    *bound[nc][0][pos[nc][0]], *bound[nc][0][pos[nc][0]-1],
                    c
                )) cnum++;
                if (clist[cnum].intersect(
                    *bound[c][e][pos[c][e]-1], *bound[c][e][pos[c][e]],
                    *bound[nc][1][pos[nc][1]], *bound[nc][1][pos[nc][1]-1],
                    c
                )) cnum++;
                // don't consider last point: because already included from other envelope
                if (pos[c][e] < boundN[c][e]-1 && clist[cnum].between(*bound[c][e][pos[c][e]],
                                        *bound[nc][1][pos[nc][1]-1],*bound[nc][1][pos[nc][1]],
                                        *bound[nc][0][pos[nc][0]-1],*bound[nc][0][pos[nc][0]],
                                        c
                                )) cnum++;
                pos[c][e]++;
            }
            break;
        } 
        else 
        {
            // determine which vertex to consider next
            if (!notline[c]) side[c]=0;
            else side[c]=((*bound[c][0][pos[c][0]])[0] < (*bound[c][1][pos[c][1]])[0])?0:1;
        }
        if ((*bound[c][side[c]][pos[c][side[c]]])[0] > (*bound[nc][side[nc]][pos[nc][side[nc]]])[0]) 
        {
            c=nc;
            nc=c?0:1;
        }
    }
    
    for (int i=0;i<cnum;i++) 
    {
        pt[index[0]]=clist[i][0];
        pt[index[1]]=clist[i][1];
        pt[index[2]]= -(clist[i][0]*n[index[0]]+ clist[i][1]*n[index[1]] - d)/n[index[2]];
        contact->addPoint(n, pt, clist[i][2], maxerr, fr);
    }

    return(cnum);
}

//---------------------------------------------------------------------------

void cDynPrimPair::findContacts(const cDynTime& time, cDynContact* contact) 
{
    static int iter=0;
    cDynVector3 localA, localB;
    cDynVector3 globalA, globalB;

    iter++;
    cDynDistPoint(primA_->v,&witness_,0,&localA);
    cDynDistPoint(primB_->v,&witness_,1,&localB);
    const cDynFrame& gA=nodeA_->globalFrame(time);
    const cDynFrame& gB=nodeB_->globalFrame(time);
    globalA.multiply(gA, localA);
    globalB.multiply(gB, localB);
    cDynVector3 ab;
    ab.subtract(globalB, globalA);
    double mag=ab.magnitude();

    double errA=primA_->err*0.5f;
    double errB=primB_->err*0.5f;
    double rA=primA_->r - errA;
    double rB=primB_->r - errB;

    double err=mag - (rA + rB);

    // find error range
    double min=(-errA > err - errB)?-errA:err - errB;
    double max=( errA < err + errB)? errA:err + errB;
    // use midpoint of range as point of contact
    double mid=0.5f*(min+max);
    double c=rA + mid;

    cDynVector3 normal;
    normal.multiply(ab, 1.0f/mag);
    cDynVector3 contactPoint;
    contactPoint.multiply(normal, c);
    contactPoint += globalA;

    //see if contacts have friction
    cDynFrictionRecord*  fr= new cDynFrictionRecord(normal,this);
    
    if (primA_->num <= 1 || primB_->num <= 1 
        || findPatch(contactPoint, normal,
         min - mid - errA, max - mid + errB, errA+errB,
         rA, rB, gA, gB, contact, fr) == 0) 
    {
        contact->addPoint(normal, contactPoint, err, errA+errB, fr);
    }
    fr->unreference();

/*
    cDynPrintf("%d: contact %s(%d) <----> %s(%d) n=[%5.9f,%5.9f,%5.9f] r=[%5.9f,%5.9f,%5.9f] err=%5.9f\n",
        iter,
        nodeA_->name()->string(),
        primA_->id,
        nodeB_->name()->string(),
        primB_->id,
        normal[0],
        normal[1],
        normal[2],
        contactPoint[0],
        contactPoint[1],
        contactPoint[2],
        err
    );
*/
}

//---------------------------------------------------------------------------


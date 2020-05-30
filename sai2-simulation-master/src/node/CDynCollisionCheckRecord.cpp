//===========================================================================
/*
    This file is part of the dynamics library.
    Copyright (C) 2014, Artificial Intelligence Laboratory,
    Stanford University. All rights reserved.

    \version   $MAJOR.$MINOR.$RELEASE $Rev: 1242 $
*/
//===========================================================================

//---------------------------------------------------------------------------
#include "node/CDynCollisionCheckRecord.h"
#include "utility/CDynLogger.h"
#include <math.h>
#include "object/CDynObject.h"
#include "node/CDynPrimPair.h"
#include "node/CDynBaseNode.h"
#include "node/CDynContact.h"
#include "distance/CDynBSphere.h"
//---------------------------------------------------------------------------
static const int CDYN_STACKSIZE=32;
//---------------------------------------------------------------------------
#define pop(a,b)  {top--;a=stack[top][0];b=stack[top][1];}
#define push(a,b) {if (top == CDYN_STACKSIZE) checkBS(nodeA,nodeB,Tab,a,b); else {stack[top][0]=a;stack[top][1]=b;top++;}}
#define pushD(a,b) {if (top == CDYN_STACKSIZE) distance(nodeA,nodeB,Tab,a,b,error,min,dis); else {stack[top][0]=a;stack[top][1]=b;top++;}}
//---------------------------------------------------------------------------

cDynPrimPair* cDynCollisionCheckRecord::distance(cDynObject* nodeA, 
    cDynObject* nodeB, 
    const cDynTransform& Tab, 
    const cDynBSphere* rootA, 
    const cDynBSphere* rootB, 
    double error, 
    double min, 
    double max)
{
    const cDynBSphere *stack[CDYN_STACKSIZE][2];
    const cDynBSphere* a;
    const cDynBSphere* b;
    cDynPrimPair* pnear=NULL;
    double dis=max;

    if (rootA == NULL) stack[0][0]=nodeA->geometry.bs();
    else stack[0][0]=rootA;
    if (rootB == NULL) stack[0][1]=nodeB->geometry.bs();
    else stack[0][1]=rootB;
    int top=1;
    cDynVector3 p;

    while (top > 0) 
    {
        pop(a,b);
        // move b to frame of a
        p.multiply(Tab, CDYN_BSPHERE_P(b));

        // test for intersection
        double radius = CDYN_BSPHERE_R(b) + CDYN_BSPHERE_R(a);
        double m = (dis + radius);
        cDynVector3 d;
        d.subtract(p, CDYN_BSPHERE_P(a));
        if (d[0] > -m && d[0] < m 
         && d[1] > -m && d[1] < m
         && d[2] > -m && d[2] < m) { 
            double dsqr=d.dot(d);
            double sdis=cDynSqrt(dsqr);
            if (sdis < m) 
            { 
                // spheres collide
                if (CDYN_BSPHERE_PRIM(a) != NULL && CDYN_BSPHERE_PRIM(b) != NULL) 
                {
                    cDynPrimPair* p=new cDynPrimPair(
                        nodeA, CDYN_BSPHERE_PRIM(a),
                        nodeB, CDYN_BSPHERE_PRIM(b)
                    );
                    if (p->checkPrim(Tab) != 2) 
                    {
                        double ndis=p->distance();
                        if (ndis <= min) return(p);
                        if (ndis < dis) 
                        {
                            dis=ndis;
                            if (pnear != NULL)
                                delete pnear;
                            pnear=p;
                        } 
                        else 
                        {
                            delete p;
                        }
                    } 
                    else 
                    {
                        delete p;
                    }
                } 
                else if (CDYN_BSPHERE_PRIM(a) != NULL) 
                {
                    biasB();
                    if (CDYN_BSPHERE_LEFT(b) != NULL) pushD(a,CDYN_BSPHERE_LEFT(b));
                    if (CDYN_BSPHERE_RIGHT(b) != NULL) pushD(a,CDYN_BSPHERE_RIGHT(b));
                } 
                else if (CDYN_BSPHERE_PRIM(b) != NULL) 
                {
                    biasA();
                    if (CDYN_BSPHERE_LEFT(a) != NULL) pushD(CDYN_BSPHERE_LEFT(a),b);
                    if (CDYN_BSPHERE_RIGHT(a) != NULL) pushD(CDYN_BSPHERE_RIGHT(a),b);
                } 
                else 
                {
                    if (CDYN_BSPHERE_R(a) > CDYN_BSPHERE_R(b)) 
                    {
                        biasA();
                        if (CDYN_BSPHERE_LEFT(a)  != NULL) pushD(CDYN_BSPHERE_LEFT(a),b);
                        if (CDYN_BSPHERE_RIGHT(a) != NULL) pushD(CDYN_BSPHERE_RIGHT(a),b);
                    } 
                    else 
                    {
                        biasB();
                        if (CDYN_BSPHERE_LEFT(b)  != NULL) pushD(a,CDYN_BSPHERE_LEFT(b));
                        if (CDYN_BSPHERE_RIGHT(b) != NULL) pushD(a,CDYN_BSPHERE_RIGHT(b));
                    }
                }
            }
        }
    }
    return(pnear);
}

//---------------------------------------------------------------------------

    int cDynCollisionCheckRecord::checkPrim(cDynObject* nodeA, 
        cDynObject* nodeB, 
        const cDynTransform& Tab, 
        const cDynPrim* a, 
        const cDynPrim* b)
{
    cDynPrimPair* p=new cDynPrimPair(
        nodeA, a,
        nodeB, b
    );

    int type=p->checkPrim(Tab);

    if (type == -1) 
    {
        penetration(p);
        if (!option(CDYN_ALL_PENETRATION))
            return(-1);
    }
    
    if (type == 0) 
    {
        contact(p);
        if (!option(CDYN_ALL_CONTACT)) 
        {
            if (np()) return(-1);
            else return(0);
        }
    }

    if (p != NULL) delete p;
    return(1);
}

//---------------------------------------------------------------------------

int cDynCollisionCheckRecord::checkBS2(cDynObject* nodeA, 
    cDynObject* nodeB, 
    const cDynTransform& Tab, 
    const cDynBSphere* rootA, 
    const cDynBSphere* rootB,
    int threshold)
{
    if (rootA == NULL) rootA=nodeA->geometry.bs();
    if (rootB == NULL) rootB=nodeB->geometry.bs();

    if (CDYN_BSPHERE_LEFT(rootA)==NULL && CDYN_BSPHERE_LEFT(rootB)==NULL) 
    {
        checkPrim(nodeA,nodeB,Tab,CDYN_BSPHERE_PRIM(rootA),CDYN_BSPHERE_PRIM(rootB));
    } 
    else 
    {
        const cDynBSphere* s[2][2];
        cDynVector3 p[2];
        cDynVector3 d;
        double radius;
        unsigned int index[4];
        double m[4];

        s[0][0]=CDYN_BSPHERE_LEFT(rootA);
        s[0][1]=CDYN_BSPHERE_RIGHT(rootA);
        s[1][0]=CDYN_BSPHERE_LEFT(rootB);
        s[1][1]=CDYN_BSPHERE_RIGHT(rootB);

        if (!s[0][0]) s[0][0]=rootA;
        if (!s[1][0]) s[1][0]=rootB;

        {for (int i=0;i<4;i++) index[i]=i;}
        {for (int i=0;i<4;i++) m[i]=0.0;}

        if (s[1][0]) p[0].multiply(Tab, CDYN_BSPHERE_P(s[1][0]));
        if (s[1][1]) p[1].multiply(Tab, CDYN_BSPHERE_P(s[1][1]));

        unsigned int n=1;
        if (s[0][1]) n += n;
        else index[1]=2;
        if (s[1][1]) n += n;

        // do test on all pairs of spheres
        {for (unsigned int i=0;i<n;i++) 
        {
            unsigned int a= index[i] & 1;
            unsigned int b= (index[i] & 2)>>1;

            // test for intersection
            radius = CDYN_BSPHERE_R(s[0][a]) + CDYN_BSPHERE_R(s[1][b]);
            d.subtract(p[b], CDYN_BSPHERE_P(s[0][a]));
            if (d[0] > -radius && d[0] < radius 
             && d[1] > -radius && d[1] < radius
             && d[2] > -radius && d[2] < radius) 
            {
                double rsqr=radius*radius;
                double dsqr=d.dot(d);
                m[i]=dsqr-rsqr;
            } 
            else 
            {
                m[i]=0.0;
            }
        }}

        // sort by amount of penetration
        // because list is so small (1,2,4 elements) i'll use bubblesort
        {for (unsigned int i=0;i < n-1; i++) 
        {
            for (unsigned int j=n-1;j >= i+1;j--) 
            {
                if (m[j] < m[j-1]) 
                { 
                    // swap
                    double tf=m[j];
                    unsigned int ti=index[j];
                    m[j]=m[j-1];index[j]=index[j-1];
                    m[j-1]=tf; index[j-1]=ti;
                }
            }
        }}
        {for (unsigned int i=0;i<n && m[i] < 0.0 && np() < threshold;i++) 
        {
            unsigned int a= index[i] & 1;
            unsigned int b= (index[i] & 2)>>1;
            checkBS2(nodeA,nodeB,Tab,s[0][a],s[1][b]);
        }}
    }

    return(np()?-1:(nc()?0:1));
}

//---------------------------------------------------------------------------

int cDynCollisionCheckRecord::checkBS(cDynObject* nodeA, 
    cDynObject* nodeB, 
    const cDynTransform& Tab, 
    const cDynBSphere* rootA, 
    const cDynBSphere* rootB)
{
    const cDynBSphere *stack[CDYN_STACKSIZE][2];
    const cDynBSphere* a;
    const cDynBSphere* b;

    if (rootA == NULL) stack[0][0]=nodeA->geometry.bs();
    else stack[0][0]=rootA;
    if (rootB == NULL) stack[0][1]=nodeB->geometry.bs();
    else stack[0][1]=rootB;
    int top=1;
    cDynVector3 p;

    while (top > 0) 
    {
        pop(a,b);
    
        // move b to frame of a
        p.multiply(Tab, CDYN_BSPHERE_P(b));

        // test for intersection
        double radius = CDYN_BSPHERE_R(b) + CDYN_BSPHERE_R(a);
        cDynVector3 d;
        d.subtract(p, CDYN_BSPHERE_P(a));
        if (d[0] > -radius && d[0] < radius 
         && d[1] > -radius && d[1] < radius
         && d[2] > -radius && d[2] < radius) 
        { 
            double rsqr=radius*radius;
            double dsqr=d.dot(d);
            if (dsqr < rsqr) 
            { 
                // spheres collide
                if (CDYN_BSPHERE_PRIM(a) != NULL && CDYN_BSPHERE_PRIM(b) != NULL) 
                {
                    cDynPrimPair* p=new cDynPrimPair(
                        nodeA, CDYN_BSPHERE_PRIM(a),
                        nodeB, CDYN_BSPHERE_PRIM(b)
                    );
/*
                    {
                        static int iter=0;
                        cDynPrintf("CHECKPRIM %d %lx, %lx, %lx, %lx\n",	iter,
                            (unsigned long)nodeA, (unsigned long)nodeB, (unsigned long)CDYN_BSPHERE_PRIM(a), (unsigned long)CDYN_BSPHERE_PRIM(b));
                        iter++;
                    }
*/
                    int type=p->checkPrim(Tab);
                    
                    if (type == -1) 
                    {
                        penetration(p);
                        if (!option(CDYN_ALL_PENETRATION))
                            return(-1);
                    }
                    
                    if (type == 0) 
                    {
                        contact(p);
                        if (!option(CDYN_ALL_CONTACT)) 
                        {
                            if (np()) return(-1);
                            else return(0);
                        }
                    }
                    if (p != NULL) delete p;

                } 
                else if (CDYN_BSPHERE_PRIM(a) != NULL) 
                {
                    biasB();
                    if (CDYN_BSPHERE_LEFT(b) != NULL) push(a,CDYN_BSPHERE_LEFT(b));
                    if (CDYN_BSPHERE_RIGHT(b) != NULL) push(a,CDYN_BSPHERE_RIGHT(b));
                } 
                else if (CDYN_BSPHERE_PRIM(b) != NULL) 
                {
                    biasA();
                    if (CDYN_BSPHERE_LEFT(a) != NULL) push(CDYN_BSPHERE_LEFT(a),b);
                    if (CDYN_BSPHERE_RIGHT(a) != NULL) push(CDYN_BSPHERE_RIGHT(a),b);
                } 
                else 
                {
                    if (CDYN_BSPHERE_R(a) > CDYN_BSPHERE_R(b)) 
                    {
                        biasA();
                        if (CDYN_BSPHERE_LEFT(a)  != NULL) push(CDYN_BSPHERE_LEFT(a),b);
                        if (CDYN_BSPHERE_RIGHT(a) != NULL) push(CDYN_BSPHERE_RIGHT(a),b);
                    } else {
                        biasB();
                        if (CDYN_BSPHERE_LEFT(b)  != NULL) push(a,CDYN_BSPHERE_LEFT(b));
                        if (CDYN_BSPHERE_RIGHT(b) != NULL) push(a,CDYN_BSPHERE_RIGHT(b));
                    }
                }
            }
        }
    }

    if (np()) return(-1);
    if (nc()) return(0);
    
    return(1);
}

//---------------------------------------------------------------------------

void cDynCollisionCheckRecord::addContacts(const cDynTime& time)
{
    cDynPrimPair* next=NULL;
    for (cDynPrimPair* p=clist_;p != NULL; p=next) 
    {
        next=p->next_;
/*
        cDynPrintf("   -- %s(%d) <--> %s(%d)\n",
            p->nodeA_->name()->string(),
            p->primA_->id,
            p->nodeB_->name()->string(),
            p->primB_->id
        );
*/
        // side effect is that the cDynPrimPair is no longer valid (thats why we use next);
        p->nodeA_->baseNode()->contact()->contact(time,p);
    }
    
    // disconnect PrimPairs from the list
    clist_=NULL;
}

//---------------------------------------------------------------------------
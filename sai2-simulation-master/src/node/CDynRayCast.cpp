//===========================================================================
/*
    This file is part of the dynamics library.
    Copyright (C) 2014, Artificial Intelligence Laboratory,
    Stanford University. All rights reserved.

    \version   $MAJOR.$MINOR.$RELEASE $Rev: 1242 $
*/
//===========================================================================

//---------------------------------------------------------------------------
#include "node/CDynRayCast.h"
#include "distance/CDynBSphere.h"
//---------------------------------------------------------------------------
#include <math.h>
//---------------------------------------------------------------------------
const int CDYN_RAYPRIM_CACHE_SIZE=13;
//---------------------------------------------------------------------------

cDynRayCast::cDynRayCast(cDynVector3& p, cDynVector3& dir, double r)
{
    p_=p;
    dir_=dir;
    dir_.normalize();
    r_=fabs(r);
    maximum_=1e99;
    reset();
    update();
}

//---------------------------------------------------------------------------

void cDynRayCast::update()
{
    int c=0;
    for (int i=1;i<3;i++) if (fabs(dir_[i]) < fabs(dir_[c])) c=i;
    cDynVector3 x,y;
    x.zero();
    x[c]=1;
    y.crossMultiply(dir_, x);
    y.normalize();
    x.crossMultiply(y, dir_);

    cDynTransform T;
    int j;
    for (j=0;j<3;j++) T.rotation(j,0)=x[j];
    for (j=0;j<3;j++) T.rotation(j,1)=y[j];
    for (j=0;j<3;j++) T.rotation(j,2)=dir_[j];
    T.translation(0) = p_[0];
    T.translation(1) = p_[1];
    T.translation(2) = p_[2];

    T_.inverse(T);
}

//---------------------------------------------------------------------------

void cDynRayCast::check(cDynObject* obj, cDynTime& time, bool tree)
{
    if (obj->geometry.bs() && !obj->isFixed()) 
    {
        cDynTransform Tnw;
        Tnw.set(obj->globalFrame(time));
        cDynTransform Tnr;
        Tnr.multiply(transform(), Tnw);

        double top;
        double d;
        double range;
        checkOne(obj->geometry.bs(), Tnr,d,range);

        if (d-range < max_) 
        {
            cDynPrim* cache[CDYN_RAYPRIM_CACHE_SIZE];
            for (int i=0;i<CDYN_RAYPRIM_CACHE_SIZE;i++) cache[i]=NULL;
            // for numerical accuracy use nearest distance known to be beyond object
            top=(((d+range)<max_)?(d+range):max_);
            double tmp=max_;max_=top;
            check(obj, obj->geometry.bs(), Tnr, cache, CDYN_RAYPRIM_CACHE_SIZE);
            if (max_ == top) max_=tmp; // did not hit anything return last good distance
        }
    }
    if (tree)
        for (cDynObject* c=obj->child();c != NULL; c=c->sibling())
            check(c, time, tree);
}

//---------------------------------------------------------------------------

void cDynRayCast::check(cDynBaseNode* base, cDynTime& time)
{
    check(base->root(), time,true);
}

//---------------------------------------------------------------------------

void cDynRayCast::check(cDynWorld* world, cDynTime& time)
{
    for (cDynBaseNode* b=world->firstBase(); b != NULL; b=b->next()) 
    {
        check(b,time);
    }
}

//---------------------------------------------------------------------------

#if 0
double cDynRayCast::checkOne(const cDynBSphere* bs, cDynTransform& Tnr)
{
    if (bs == NULL) return(max_);
    cDynVector3 p;
    p.multiply(Tnr, CDYN_BSPHERE_P(bs));

    if (p[2] - CDYN_BSPHERE_R(bs) >= max_) return(maximum_); // to far
    if (p[2] + CDYN_BSPHERE_R(bs) <= min_) return(maximum_); // to close
    double r2=CDYN_BSPHERE_R(bs) + r_;
    if (p[0] > r2 || p[0] < -r2) return(maximum_); // to the side
    if (p[1] > r2 || p[1] < -r2) return(maximum_); // to the side
    double dsqr=p[0]*p[0]+p[1]*p[1];

    if (dsqr > r2*r2) return(maximum_);
    return(p[2] + CDYN_BSPHERE_R(bs));
}
#endif

void cDynRayCast::checkOne(const cDynBSphere* bs, cDynTransform& Tnr, double& z, double& range)
{
    if (bs == NULL) { z=maximum_;range=0.0; return; }
    cDynVector3 p;
    p.multiply(Tnr, CDYN_BSPHERE_P(bs));

    double r2=CDYN_BSPHERE_R(bs) + r_;
    if (p[2] - CDYN_BSPHERE_R(bs) >= max_  // to far
      || p[2] + CDYN_BSPHERE_R(bs) <= min_ // to close
      || p[0] > r2 || p[0] < -r2        // to the side
      || p[1] > r2 || p[1] < -r2        // to the other side
    ) { z=maximum_;range=0.0; return; }
    double dsqr=p[0]*p[0]+p[1]*p[1];

    if (dsqr > r2*r2) { z=maximum_;range=0.0; return; }
    z=p[2];
    range=sqrt(r2*r2 - dsqr);
    return;
}

//---------------------------------------------------------------------------

void cDynRayCast::check(cDynObject* obj, const cDynBSphere* bs, cDynTransform& Tnr, cDynPrim* cache[], int n)
{
    if (bs == NULL) return;
    if (CDYN_BSPHERE_LEFT(bs) == NULL && CDYN_BSPHERE_RIGHT(bs) == NULL) 
    {
        int h=CDYN_BSPHERE_PRIM(bs)->id%n;
        if (cache[h] == CDYN_BSPHERE_PRIM(bs)) return;
        check(obj, CDYN_BSPHERE_PRIM(bs), Tnr);
        cache[h]=CDYN_BSPHERE_PRIM(bs);
    } 
    else 
    {
        double z,range;
        checkOne(CDYN_BSPHERE_LEFT(bs),Tnr,z,range);
        double dl=z-range;
        checkOne(CDYN_BSPHERE_RIGHT(bs),Tnr,z,range);
        double dr=z-range;
        if (dl < dr) 
        {
            if (dl < max_) check(obj,CDYN_BSPHERE_LEFT(bs),Tnr,cache,n);
            if (dr < max_) check(obj,CDYN_BSPHERE_RIGHT(bs),Tnr,cache,n);
        } 
        else 
        {
            if (dr < max_) check(obj,CDYN_BSPHERE_RIGHT(bs),Tnr,cache,n);
            if (dl < max_) check(obj,CDYN_BSPHERE_LEFT(bs),Tnr,cache,n);
        }
    }
}

//---------------------------------------------------------------------------

void cDynRayCast::check(cDynObject* obj, cDynPrim* prim, cDynTransform& Tnr)
{
    cDynPrim ray=cDynPrim(2);
    cDynTransform Tinv;
    Tinv.inverse(Tnr);

    cDynVector3 zl;
    zl.set(0,0,min_);
    cDynVector3 zh;
    zh.set(0,0,max_);
    ray.v[0].multiply(Tinv, zl);
    ray.v[1].multiply(Tinv, zh);

    cDynDistWitness witness;
    CDYN_DISTW_NOV(&witness)=0;
    double dis=cDynDistDistance(prim->v,prim->num,ray.v,ray.num,&witness,0,false);
    //double dis=cDynDistDistance2(prim->v,prim->num,ray.v,ray.num,&witness,0,false);

    dis -= prim->r + r_;

    if (dis > 0.0f) return;
    if (CDYN_DISTW_NOV(&witness) < 2) return;
    cDynVector3 p;
    cDynDistPoint(prim->v,&witness,0,&p);
    cDynVector3 pr;
    pr.multiply(Tnr, p);
    double dr=cDynSqrt(pr[0]*pr[0] + pr[1]*pr[1]);
    double di=(dr >= prim->r)?pr[2]:pr[2]-cDynSqrt(prim->r*prim->r - dr*dr);

    if (di < max_) 
    {
        if (di < min_) di= min_;
        node_=obj;
        prim_=prim;
        max_=di;
    }
}

//---------------------------------------------------------------------------

//===========================================================================
/*
    This file is part of the dynamics library.
    Copyright (C) 2014, Artificial Intelligence Laboratory,
    Stanford University. All rights reserved.

    \version   $MAJOR.$MINOR.$RELEASE $Rev: 1242 $
*/
//===========================================================================

//---------------------------------------------------------------------------
#include <math.h>
#include "object/CDynObject.h"
#include "node/CDynContactPointCache.h"
//---------------------------------------------------------------------------
cDynContactPointCache g_dynTheContactPointCache;
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
/*!
    Compute a hash function to place a given contact point randomly in 
    the cache.
*/
//---------------------------------------------------------------------------
unsigned int cDynContactPointCache::h(const cDynContact* c, 
    const cDynVector3& n, 
    const cDynVector3& r, 
    const double err, 
    const double maxerr)
{
    const unsigned int sn[3]={5,7,13};
    const unsigned int sr[3]={3,23,29};
    const unsigned int se=31;

    cDynObject* a=c->prim_->nodeA_;
    cDynObject* b=c->prim_->nodeB_;
    unsigned int h= (a->uid() + b->uid())%CDYN_POINTCACHESIZE;
    unsigned int hn = 0;
    unsigned int j = 0;

    for (int i=0;i<3;i++) 
    {
        j=(unsigned int)((0.5f*(n[i]+1.0f))*CDYN_POINTCACHESIZE);
        hn +=sn[i]*j;
    }

    for (int i=0;i<3;i++) 
    {
        j=(unsigned int)(fabs((1e4*r[i])+0.5f));
        hn +=sr[i]*j;
    }

    j = (unsigned int)(0.5f*(err + maxerr)*CDYN_POINTCACHESIZE);
    hn +=se*j;
    h = (h + hn)%CDYN_POINTCACHESIZE;
    return(h);
}

//---------------------------------------------------------------------------

int cDynContactPointCache::check(const cDynContact* c, 
    const cDynVector3& n,
    const cDynVector3& r, 
    const double err, 
    const double maxerr)
{
    unsigned int hash=h(c,n,r,err,maxerr);
    if (num_[hash] != callnum_) {miss_++; return(hash);};
    cDynContactPoint* p=point_[hash];
    if (c->prim_->nodeA_ != p->contact_->prim_->nodeA_) return(hash);
    if (c->prim_->nodeB_ != p->contact_->prim_->nodeB_) return(hash);
    if (fabs(err - p->err_) > 1e-3*maxerr) return(hash);
    cDynVector3 e;
    e.subtract(r, p->r_);
    if (e.magnitude() > 1e-3*maxerr) return(hash);
    e.subtract(n, p->s_[0]);
    if (e.magnitude() > 1e-6) return(hash);

    // looks like we have a hit
    hit_++;
    return(-1);
}

//---------------------------------------------------------------------------

void cDynContactPointCache::insert(cDynContactPoint* p, const int h)
{
    if (num_[h] == callnum_) collision_++;
    num_[h]=callnum_;
    point_[h]=p;
}

//---------------------------------------------------------------------------
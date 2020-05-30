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
#include "distance/CDynBSphere.h"
#include "distance/CDynPrim.h"
//---------------------------------------------------------------------------
#define CDYN_CACHESIZE 7
//---------------------------------------------------------------------------

cDynBSphere::cDynBSphere(cDynVector3* v, double radius, cDynPrim *p)
{
    CDYN_BSPHERE_P(this) = *v;
    CDYN_BSPHERE_R(this) = radius;
    CDYN_BSPHERE_PRIM(this) = p;
    CDYN_BSPHERE_STATE(this) = 0;
    CDYN_BSPHERE_LEFT(this) = NULL;
    CDYN_BSPHERE_RIGHT(this) = NULL;
}

//---------------------------------------------------------------------------

void cDynBSphere::remove()
{
    if (CDYN_BSPHERE_LEFT(this) != NULL) CDYN_BSPHERE_LEFT(this)->remove();
    if (CDYN_BSPHERE_RIGHT(this) != NULL) CDYN_BSPHERE_RIGHT(this)->remove();
    delete this;
}

//---------------------------------------------------------------------------

double cDynBSphere::max(const cDynVector3 *dir, 
    cDynPrim *wprim, 
    int *witness, 
    const double d, 
    cDynPrim* cache[]) const
{
    const cDynBSphere *bs = this;
    double nd=d;
    if (cache == NULL) 
    {
        cDynPrim* c[CDYN_CACHESIZE];
        for (int i=0;i<CDYN_CACHESIZE;i++) c[i]=NULL;
        cache = c;
        if (wprim != NULL) nd=dir->dot(wprim->v[*witness])+wprim->r;
        else nd=dir->dot(CDYN_BSPHERE_P(bs))-CDYN_BSPHERE_R(bs);
    }
    if (CDYN_BSPHERE_PRIM(bs) != NULL) 
    {
        int h=CDYN_BSPHERE_PRIM(bs)->id%CDYN_CACHESIZE;
        if (cache[h] != CDYN_BSPHERE_PRIM(bs)) 
        {
            cache[h]=CDYN_BSPHERE_PRIM(bs);
            for (int i=0;i<CDYN_BSPHERE_PRIM(bs)->num;i++) 
            {
                double m=dir->dot(CDYN_BSPHERE_PRIM(bs)->v[i]) + CDYN_BSPHERE_PRIM(bs)->r;
                if (m > nd) 
                {
                    nd=m;
                    *witness = i;
                    wprim=CDYN_BSPHERE_PRIM(bs);
                }
            }
        }
    } 
    else 
    {
        cDynBSphere* b[2]={CDYN_BSPHERE_LEFT(bs),CDYN_BSPHERE_RIGHT(bs)};
        double m[2]={d,d};
        int i;
        for (i=0;i<2;i++) 
        {
            if (b[i] != NULL) 
            {
                m[i]=dir->dot(CDYN_BSPHERE_P(b[i])) + CDYN_BSPHERE_R(b[i]);
            }
        }
        if (m[1] > m[0]) 
        {
            cDynBSphere* bt=b[0]; b[0]=b[1]; b[1]=bt;
            double    mt=m[0]; m[0]=m[1]; m[1]=mt;
        }
        for (i=0;i<2;i++) 
        {
            if (b[i] != NULL && m[i] > nd) 
            {
                nd=b[i]->max(dir,wprim,witness,nd,cache);
            }
        }
    }

    return(nd);
}

//---------------------------------------------------------------------------
//===========================================================================
/*
    This file is part of the dynamics library.
    Copyright (C) 2014, Artificial Intelligence Laboratory,
    Stanford University. All rights reserved.

    \version   $MAJOR.$MINOR.$RELEASE $Rev: 1242 $
*/
//===========================================================================

//---------------------------------------------------------------------------
#include <assert.h>
#include "distance/CDynBSphere.h"
#include "distance/CDynPrim.h"
#include "distance/CDynCover.h"
#include "object/CDynObject.h"
//---------------------------------------------------------------------------
//#define CDYN_RADIUS 0.01f;
#define CDYN_RADIUS 0.09f;
//#define CDYN_RADIUS 2.0f;
//---------------------------------------------------------------------------
static double _cDynRadius = CDYN_RADIUS;
//---------------------------------------------------------------------------

static void _cDynCoverLine(cDynGeometry *geom, 
    cDynVector3& p, 
    cDynVector3& q, 
    cDynPrim *prim, 
    double rmin, 
    double width, 
    double rmax)
{
    double w2=2.0f*width;
    cDynVector3 dp;
    dp.subtract(q, p);
    double d=dp.magnitude();

    if (width > rmin) 
    {
        double dm=cDynSqrt(d*d*0.25f+width*width);
        if ((rmin + dm) < rmax) 
        {
            cDynVector3 c;
            c.add(p,q);
            c*=0.5f;
            geom->bsphere(new cDynBSphere(&c,dm + rmin,prim));
            return;
        }
    } 
    else 
    {
        if ((rmin + d + rmin) < rmax + rmax) 
        {
            cDynVector3 c;
            c.add(p,q);
            c*=0.5f;
            geom->bsphere(new cDynBSphere(&c,d*0.5f + rmin,prim));
            return;
        } 
        else 
        {
            width=rmax/(double)M_SQRT2;
        }
    }
    cDynVector3 dn;
    dn.multiply(dp, 1.0f/d);
    double s=2.0f*(width - rmin);
    cDynVector3 tmpV;
    tmpV.multiply(dn, s);
    dp -= tmpV;
    d -= s;

    int n = (int)(d/w2);
    if (n*w2 < d || n==0 ) n++;

    cDynVector3 step;
    step.multiply(dp, 1.0f/n);
    cDynVector3 c;
    c.multiply(dn, width - rmin);
    c+=p;

    for (int i=0;i<=n;i++) 
    {
        geom->bsphere(new cDynBSphere(&c,rmax,prim));
        c += step;
    }
}

//---------------------------------------------------------------------------

void cDynCoverRadius(const double r)
{
    assert(r > 0.0f);
    _cDynRadius=r;
}

//---------------------------------------------------------------------------

void cDynCoverHull(cDynGeometry *geom, cDynPrim *prim)
{
    // for now if it is not very simple geometry surround it with just one sphere
    switch (prim->num) 
    {
        case 0: return;
        case 1: geom->bsphere(new cDynBSphere(&prim->v[0],prim->r,prim)); return;
        case 2: _cDynCoverLine(geom,prim->v[0], prim->v[1],prim,prim->r,prim->r,prim->r*(double)M_SQRT2); return;
    }
    int num=prim->num;

    // find centroid
    cDynVector3 avg=prim->v[0];

    for (int i=1;i<num;i++) 
    {
        avg += prim->v[i];
    }
    avg *= (1.0f/(double)num);

    // find Enclosing Sphere
    double rsqr=0.0f;
    cDynVector3 d;
    for (int i=0;i<num;i++) 
    {
        d.subtract(prim->v[i], avg);
        double rr=d.dot(d);
        if (rr > rsqr) rsqr=rr;
    }
    
    double r=cDynSqrt(rsqr) + prim->r;

    geom->bsphere(new cDynBSphere(&avg,r,prim));

    return;
}

//---------------------------------------------------------------------------

void cDynCoverPoly(cDynGeometry *geom, cDynPrim *prim)
{
    double rs;

    if (prim->r < _cDynRadius) rs=_cDynRadius;
    else rs=prim->r;
    switch (prim->num) 
    {
        case 0: return;
        case 1: geom->bsphere(new cDynBSphere(&prim->v[0],prim->r,prim)); return;
        case 2: _cDynCoverLine(geom,prim->v[0], prim->v[1],prim,prim->r,prim->r,rs*(double)M_SQRT2); return;
    }
    int num=prim->num;
    double rmin=prim->r;
    double width=rs;
    double rmax=rs*(double)M_SQRT2;

    double w2=2.0f*width;

    // find centroid to see if poly will fit in one sphere
    // and determine how far polygon deviates from defining plane

    double err,maxerr=0;
    cDynVector3& n=prim->normal;
    double nd=prim->d;
    cDynVector3 avg,max,min;
    avg=prim->v[0];
    max=avg;
    min=max;

    for (int i=1;i<num;i++) 
    {
        avg += prim->v[i];
        err = cDynFabs(n.dot(prim->v[i])+nd);
        min.minimum(prim->v[i]);
        max.maximum(prim->v[i]);
        if (err > maxerr) maxerr=err;
    }
    avg *= (1.0f/(double)num);
    if (maxerr > width) 
    {
        width = maxerr;
        rmax = width*(double)M_SQRT2;
    }

    // Find Enclosing Sphere
    double rsqr=0.0f;
    for (int i=0;i<num;i++) 
    {
        cDynVector3 d;
        d.subtract(prim->v[i], avg);
        double rr=d.dot(d);
        if (rr > rsqr) rsqr=rr;
    }
    
    double r=cDynSqrt(rsqr) + rmin;

    if (r < rmax) 
    { 
        // polygon fits inside sphere
        geom->bsphere(new cDynBSphere(&avg,r,prim)); return;
    }

    // find projection axis z for scan conversion
    cDynVector3 range;
    range.subtract(max, min);
    int x,y,z;
    if (cDynFabs(n[0]) > cDynFabs(n[1]) && cDynFabs(n[0]) > cDynFabs(n[2])) 
    {
        if (range[1] > range[2]) { y=1;x=2; }
        else { x=1; y=2; }
        z=0;
    } 
    else if (cDynFabs(n[1]) > cDynFabs(n[2])) 
    {
        if (range[0] > range[2]) { y=0;x=2; }
        else { x=0; y=2; }
        z=1;
    } 
    else 
    {
        if (range[0] > range[1]) { y=0;x=1; }
        else { x=0; y=1; }
        z=2;
    }

    // find vertex that is largest in y direction
    double s=prim->v[0][y];
    int p0,p1,i0=0,i1;
    for (int i=1;i<num;i++) 
    {
        if (prim->v[i][y] > s) 
        {
            s = prim->v[i][y];
            i0=i;
        }
    }
    s += rmin;
    p0=p1=i1=i0;

    double d= range[y] + 2.0f*rmin*cDynFabs(n[z]);
    int rows = (int)(d/w2);
    if (rows*w2 < d || rows == 0) rows++;
    double ystep = d/rows;
    double smid = s - 0.5f*ystep;

    s -= ystep;

    double xmin,xmax;
    xmin=xmax=prim->v[p0][x];

    int top=1;
    for (int i=0;i<rows;i++) 
    {
        while (prim->v[i0][y] > s && ((i0 != i1) || top)) 
        {
            p0=i0;
            if (--i0 < 0)  i0=num-1;
            top=0;
            if (prim->v[p0][x] > xmax) xmax=prim->v[p0][x];
            if (prim->v[p0][x] < xmin) xmin=prim->v[p0][x];
        }
        while (prim->v[i1][y] > s && ((i0 != i1) || top)) 
        {
            p1=i1;
            if (++i1 >= num) i1=0;
            top=0;
            if (prim->v[p1][x] > xmax) xmax=prim->v[p1][x];
            if (prim->v[p1][x] < xmin) xmin=prim->v[p1][x];
        }

        double l,xint0,xint1;

        l=prim->v[i0][y] - prim->v[p0][y];
        if (s < prim->v[i0][y] || cDynFabs(l) < 1e-8f)
            xint0=prim->v[i0][x];
        else
            xint0=xint0=prim->v[p0][x] + (prim->v[i0][x] - prim->v[p0][x])*(s-prim->v[p0][y])/l;

        l=prim->v[i1][y] - prim->v[p1][y];
        if (s < prim->v[i1][y] || cDynFabs(l) < 1e-8f)
            xint1=prim->v[i1][x];
        else
            xint1=prim->v[p1][x] + (prim->v[i1][x] - prim->v[p1][x])*(s-prim->v[p1][y])/l;

        if (xint0 > xmax) xmax=xint0;
        if (xint0 < xmin) xmin=xint0;
        if (xint1 > xmax) xmax=xint1;
        if (xint1 < xmin) xmin=xint1;

        cDynVector3 a,b;
        a[x]=xmin; a[y]=smid; a[z]=-(n[x]*xmin+n[y]*smid + nd)/n[z];
        b[x]=xmax; b[y]=smid; b[z]=-(n[x]*xmax+n[y]*smid + nd)/n[z];

        _cDynCoverLine(geom, a,b, prim, rmin, width, rmax);

        xmin=xmax=xint0;

        if (xint1 > xmax) xmax=xint1;
        if (xint1 < xmin) xmin=xint1;

        s -= ystep;
        smid -= ystep;
    }
}

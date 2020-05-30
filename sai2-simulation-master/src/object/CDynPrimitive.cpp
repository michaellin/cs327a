//===========================================================================
/*
    This file is part of the dynamics library.
    Copyright (C) 2014, Artificial Intelligence Laboratory,
    Stanford University. All rights reserved.

    \version   $MAJOR.$MINOR.$RELEASE $Rev: 1242 $
*/
//===========================================================================

//---------------------------------------------------------------------------
#include <cmath>
#include <string.h> // memcpy
#include "distance/CDynPrim.h"
#include "object/CDynObject.h"
#include "utility/CDynError.h"
#include "distance/CDynCover.h"
#include "utility/CDynLogger.h"
//---------------------------------------------------------------------------
//unsigned int cDynPrim::lastid = 0;
//---------------------------------------------------------------------------

cDynPrimitive::cDynPrimitive(cDynGeometry* geom, const cDynPrimitiveType type)
    : geom_(geom), type_(type), radius_(0.0f), error_(DEPRIMITIVE_MIN_ERROR), pvertex_(NULL), 
    nv_(0), nvallocated_(0), nchunk_(12)
{
    material_ = geom->material();
    data_=NULL;
}

//---------------------------------------------------------------------------

cDynPrimitive::~cDynPrimitive()
{
    if (pvertex_)
    {
        delete [] pvertex_;
    }
}

//---------------------------------------------------------------------------

/*
void cDynPrimitive::vertex(const cDynVector3& v)
{
    if (nv_ >= nvallocated_) 
    {
        cDynVector3 *tmp = NULL;
        int i;

        // XXX why is the size of cDynVector3 24 bytes?
        if (pvertex_)
        {
            tmp = new cDynVector3[nv_];
            for (i = 0; i < nv_; i++)
                tmp[i] = pvertex_[i];
            delete [] pvertex_;
        }

        nvallocated_ += nchunk_;
        pvertex_ = new cDynVector3[nvallocated_];

        if (tmp)
        {
            for (i = 0; i < nv_; i++)
                pvertex_[i] = tmp[i];
        }
    }

    if (geom_->frame() == NULL)
        pvertex_[nv_++]=v;
    else
        pvertex_[nv_++].multiply(geom_->frame()->top(), v);
}
*/

//---------------------------------------------------------------------------

void cDynPrimitive::vertex(const cDynVector3& v)
{
    if (nv_ >= nvallocated_) 
    {
        cDynVector3 *tmp = NULL;

        // XXX why is the size of cDynVector3 24 bytes?
        if (pvertex_)
        {
            tmp = pvertex_;
        }

        nvallocated_ += nchunk_;
        pvertex_ = new cDynVector3[nvallocated_];

        if (tmp)
        {
            memcpy(pvertex_, tmp, sizeof(cDynVector3) * nv_);
            delete [] tmp;
        }
    }

    if (geom_->frame() == NULL)
        pvertex_[nv_++]=v;
    else
        pvertex_[nv_++].multiply(geom_->frame()->top(), v);
}

//---------------------------------------------------------------------------

void cDynPrimitive::vertex(const double x, const double y, const double z = 0.0f)
{ 
    cDynVector3 vec;
    vec.set(x, y, z);
    vertex(vec);
}

//---------------------------------------------------------------------------

void cDynPrimitive::vertex(const double v[3])
{ 
    cDynVector3 vec;
    vec.set(v[0], v[1], v[2]);
    vertex(vec);
}

//---------------------------------------------------------------------------

void cDynPrimitive::radius(double r)
{
    if (r >= 0.0f)
        radius_ = r;
    else 
        cDynError(CDYN_ERROR_PRIMITIVE_NEG_RADIUS);
}

//---------------------------------------------------------------------------

void cDynPrimitive::error(double err)
{
    if (err >= DEPRIMITIVE_MIN_ERROR)
        error_ = err;
    else 
        cDynError(CDYN_ERROR_PRIMITIVE_ERROR_TOO_SMALL);
}

//---------------------------------------------------------------------------

void cDynPrimitive::point(int i)
{
    // allocate prim to hold one point
    cDynPrim *p = new cDynPrim(1);
    p->data(data());

    p->v[0] = pvertex_[i];

    // normal not defined for points 
    p->normal.zero();
    p->d = 0.0f;

    // set radius and error
    p->r = radius_ + error_;
    p->err = 2.0f * error_;

    p->m = material();

    geom_->prim(p);

    // cover the point
    cDynCoverPoly(geom_, p);
}

//---------------------------------------------------------------------------

void cDynPrimitive::line(int i0, int i1)
{
    // allocate prim to hold two points
    cDynPrim *p = new cDynPrim(2);
    p->data(data());

    p->v[0] = pvertex_[i0];
    p->v[1] = pvertex_[i1];

    // surface plane not defined for lines
    p->normal.zero();
    p->d = 0.0f;

    // set radius and error
    p->r = radius_ + error_;
    p->err = 2.0f *error_;

    p->m = material();

    geom_->prim(p);

    // cover the line
    cDynCoverPoly(geom_, p);
}

//---------------------------------------------------------------------------

void cDynPrimitive::triangle(int i0, int i1, int i2)
{
    cDynPrim *p = new cDynPrim(3);
    p->data(data());

    p->v[0] = pvertex_[i0];
    p->v[1] = pvertex_[i1];
    p->v[2] = pvertex_[i2];

    cDynVector3 v10;
    v10.subtract(p->v[1], p->v[0]);
    cDynVector3 v20;
    v20.subtract(p->v[2], p->v[0]);
    p->normal.crossMultiply(v10, v20);
    p->normal.normalize();
    p->d= -(p->normal.dot(p->v[0]));

    // set radius and error
    p->r = radius_ + error_;
    p->err = 2.0f  *error_;

    p->m = material();

    geom_->prim(p);

    // cover the line
    cDynCoverPoly(geom_, p);
}

//---------------------------------------------------------------------------

void cDynPrimitive::quad(int i0, int i1, int i2, int i3)
{
    int n = 4;
    cDynPrim *p = new cDynPrim(n);
    p->data(data());
    
    p->v[0] = pvertex_[i0];
    p->v[1] = pvertex_[i1];
    p->v[2] = pvertex_[i2];
    p->v[3] = pvertex_[i3];

    // calculate the surface normal
    double a, b, c;
    cDynVector3 mean;
    
    a = b = c = 0.0f;
    mean.zero();

    for (int i = 0; i < n; i++) 
    {
        int j;

        if (i == n - 1) 
            j = 0;
        else 
            j = i + 1;

        a += (p->v[i][1] - p->v[j][1]) * (p->v[i][2] + p->v[j][2]);
        b += (p->v[i][2] - p->v[j][2]) * (p->v[i][0] + p->v[j][0]);
        c += (p->v[i][0] - p->v[j][0]) * (p->v[i][1] + p->v[j][1]);
        mean += p->v[i];
    }

    mean *= 1.0f / n;
    p->normal.set(a, b, c);
    double d = p->normal.magnitude();
    if (d < 1e-8) 
    { 
        // Degenerate polygon need to do something FIX XXX
        delete p;
        return;
    }
    p->normal *= 1.0f / d;
    
    p->d= -(p->normal.dot(mean));

    // set radius and error
    p->r= radius_ + error_;
    p->err= 2.0f * error_;

    p->m = material();

    geom_->prim(p);

    // cover the line
    cDynCoverPoly(geom_, p);
}

//---------------------------------------------------------------------------

void cDynPrimitive::polygon(int istart, int iend)
{
    int i, n = iend - istart;
    
    cDynPrim *p = new cDynPrim(n);
    p->data(data());
    
    for (i = 0; i < n; i++) 
        p->v[i] = pvertex_[istart + i];

    // calculate the surface normal
    double a, b, c;
    cDynVector3 mean;
    
    a = b = c = 0.0f;
    mean.zero();

    for (i = 0; i < n; i++) 
    {
        int j;

        if (i == n - 1)
            j = 0;
        else 
            j = i + 1;

        a += (p->v[i][1] - p->v[j][1]) * (p->v[i][2] + p->v[j][2]);
        b += (p->v[i][2] - p->v[j][2]) * (p->v[i][0] + p->v[j][0]);
        c += (p->v[i][0] - p->v[j][0]) * (p->v[i][1] + p->v[j][1]);
        mean += p->v[i];
    }

    mean *= 1.0f / n;
    p->normal.set(a, b, c);
    double d = p->normal.magnitude();

    if (d < 1e-8) 
    { 
        // Degenerate polygon need to do something FIX XXX
        delete p;
        return;
    }
    p->normal *= 1.0f / d;
    
    p->d = -(p->normal.dot(mean));

    // set radius and error
    p->r = radius_ + error_;
    p->err = 2.0f * error_;

    p->m = material();

    geom_->prim(p);

    // cover the polygon
    cDynCoverPoly(geom_, p);
}

//---------------------------------------------------------------------------

void cDynPrimitive::hull(int istart, int iend)
{
    int i, n = iend - istart;
    
    cDynPrim *p = new cDynPrim(n);
    p->data(data());
    
    for (i = 0; i < n; i++) 
        p->v[i] = pvertex_[istart + i];

    p->normal.zero();
    p->d = 0.0;

    // set radius and error
    p->r = radius_ + error_;
    p->err = 2.0f * error_;

    p->m = material();

    geom_->prim(p);

    // cover the hull
    cDynCoverHull(geom_, p);
}

//---------------------------------------------------------------------------

void cDynPrimitive::bound(int istart, int iend)
{
    int i, n = iend - istart;
    
    cDynPrim *p = new cDynPrim(n);
    p->data(data());
    
    for (i = 0; i < n; i++) 
        p->v[i] = pvertex_[istart + i];

    p->normal.set(0, 0, 0);
    p->d = 0;

    // set radius and error
    p->r = radius_ + error_;
    p->err = 2 * error_;

    p->m = material();

    geom_->bound(p);
}

//---------------------------------------------------------------------------

void cDynPrimitive::convert(void)
{
    int i;

    if (nv_ < 1) 
        return;

    switch (type_) 
    {
    case CDYN_POINTS:
        for (i = 0; i < nv_; i++)
            point(i);
        break;

    case CDYN_LINES:
        nv_ -= nv_%2; // make sure we have an even number of points
        for (i = 0; i < nv_; i += 2)
            line(i,i+1);
        break;

    case CDYN_LINE_STRIP:
        for (i = 0; i < nv_ - 1; i++)
            line(i,i+1);
        break;

    case CDYN_LINE_LOOP:
        for (i = 0; i < nv_ - 1; i++)
            line(i,i+1);
        if (nv_ > 2) 
            line(nv_ - 1, 0);
        break;

    case CDYN_TRIANGLES:
        nv_ -= nv_ % 3; // num points must be multiple of 3
        if (nv_ == 0) 
            break;
        for (i = 0; i < nv_; i += 3)
            triangle(i, i + 1, i + 2);
        break;

    case CDYN_TRIANGLE_STRIP:
        for (i = 0; i < nv_ - 2; i++)
            triangle(i,i+1,i+2);
        break;

    case CDYN_TRIANGLE_FAN:
        for (i = 1; i < nv_ - 1; i++)
            triangle(0,i,i+1);
        break;

    case CDYN_QUADS:
        nv_ -= nv_ % 4; // num points must be multiple of 4
        if (nv_ == 0) 
            break;
        for (i = 0; i <= nv_ - 4; i += 4)
            quad(i, i + 1, i + 2, i + 3);
        break;

    case CDYN_QUAD_STRIP:
        nv_ -= nv_ % 2;// num points must be multiple of 2
        for (i = 0; i <= nv_ - 4 ; i += 2)
            quad(i, i + 1, i + 3, i + 2);
        break;

    case CDYN_POLYGON:
        switch (nv_) 
        {
        case 1:	point(0); break;
        case 2:	line(0, 1); break;
        default: polygon(0, nv_);
        }
        break;

    case CDYN_POLYGON_CONCAVE:
        break;

    case CDYN_HULL:
        hull(0,nv_);
        break;

    case CDYN_BOUND:
        bound(0,nv_);

    default:
        break;
    }
}

//---------------------------------------------------------------------------

void cDynPrimitive::end(void)
{
    geom_->end(this);
}

//---------------------------------------------------------------------------
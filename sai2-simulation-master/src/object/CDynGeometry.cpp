//===========================================================================
/*
    This file is part of the dynamics library.
    Copyright (C) 2014, Artificial Intelligence Laboratory,
    Stanford University. All rights reserved.

    \version   $MAJOR.$MINOR.$RELEASE $Rev: 1242 $
*/
//===========================================================================

//---------------------------------------------------------------------------
#include "object/CDynObject.h"
#include "distance/CDynHierarchy.h"
#include "distance/CDynBSphere.h"
#include "distance/CDynPrim.h"
//---------------------------------------------------------------------------
#include <assert.h>
//---------------------------------------------------------------------------

void cDynGeometryNode::remove()
{
    if (bs_ != NULL) bs_->remove();
    if (prim_ != NULL) prim_->remove();
    if (bound_ != NULL) bound_->remove();
}

//---------------------------------------------------------------------------

void cDynGeometryNode::build()
{
    // set new geometry
    bs_=cDynHierarchy(bs_);

    // tighten the top sphere as much as possible
    if (bs_ != NULL) 
    {
        double maxr=0.0f;
        cDynVector3& c=CDYN_BSPHERE_P(bs_);
        cDynVector3  v;
        for (cDynPrim* p=prim_;p!=NULL;p=p->next) 
        {
            for (int i=0;i<p->num;i++) 
            {
                v.subtract(c, p->v[i]);
                double r=v.magnitude() + p->r;
                if (r > maxr) maxr=r;
            }
        }
        assert(maxr < CDYN_BSPHERE_R(bs_)+1e-8);
        
        //cDynPrintf("Savings %0.5g\n",maxr/CDYN_BSPHERE_R(bs_));
        CDYN_BSPHERE_R(bs_)=maxr;
    }
}

//---------------------------------------------------------------------------

void cDynGeometryNode::prim(cDynPrim *p)
{
    // add prim to linked list
    p->next=prim_;
    prim_=p;
}

//---------------------------------------------------------------------------

void cDynGeometryNode::bound(cDynPrim *p)
{
    if (bound_ != NULL) free(bound_);
    bound_=p;
}

//---------------------------------------------------------------------------

void cDynGeometryNode::bsphere(cDynBSphere *bs)
{
    CDYN_BSPHERE_LEFT(bs)=bs_;
    CDYN_BSPHERE_RIGHT(bs)=NULL;

    bs_=bs;
}

//---------------------------------------------------------------------------

cDynGeometry::cDynGeometry(void)
{
    fs_=NULL;

    // pointers to referenced geometry
    ptr_=new cDynGeometryNode();

    // default material
    material_=NULL;

    // on creation bounding box undefined
    box_=NULL;
}

//---------------------------------------------------------------------------

cDynGeometry::cDynGeometry(cDynGeometry& src)
{
    fs_=NULL;

    material_=src.material_;
    box_=NULL;

    src.ptr_->reference();
    ptr_=src.ptr_;
}

//---------------------------------------------------------------------------

cDynGeometry& cDynGeometry::operator=(cDynGeometry& src)
{
    // reference first in case its a this = this situation
    src.ptr_->reference();
    ptr_->unreference();
    ptr_=src.ptr_;
    return *this;
}

//---------------------------------------------------------------------------

cDynGeometry::~cDynGeometry()
{
    // delete any geometry
    ptr_->unreference();
}

//---------------------------------------------------------------------------

cDynPrimitive* cDynGeometry::begin(const cDynPrimitiveType type)
{
    cDynPrimitive* p=new cDynPrimitive(this,type);
    
    return(p);
}

//---------------------------------------------------------------------------

void cDynGeometry::end(cDynPrimitive* p)
{
    // convert to prim and cover object
    p->convert();

    // free up memory
    delete p;
}

//---------------------------------------------------------------------------

void cDynGeometry::end(void)
{
    // build hierarchy
    ptr_->build();
}

//---------------------------------------------------------------------------
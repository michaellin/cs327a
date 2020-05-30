//===========================================================================
/*
    This file is part of the dynamics library.
    Copyright (C) 2014, Artificial Intelligence Laboratory,
    Stanford University. All rights reserved.

    \version   $MAJOR.$MINOR.$RELEASE $Rev: 1242 $
*/
//===========================================================================

//---------------------------------------------------------------------------
#include "matrix/CDynVector3.h"
#include "object/CDynObject.h"
#include "node/CDynContact.h"
#include "node/CDynContactPointCache.h"
#include "node/CDynConstraint.h"
//---------------------------------------------------------------------------
cDynContact* cDynContact::free_=NULL;
//---------------------------------------------------------------------------

void cDynContact::addPoint(const cDynVector3& n,
    const cDynVector3& r, 
    const double err, 
    const double maxerr, 
    cDynFrictionRecord* fr)
{
    assert(type_ == CDYN_CONTACT);
    // check if point is already in the cache
    int h=g_dynTheContactPointCache.check(this,n,r,err,maxerr);
    if (h == -1) return; // point redundant

    cDynContactPoint* cp=NULL;
    if (!fr->friction()) 
    { 
        // add constraint normal to surface
        cp=new cDynContactPoint(this,n,r,err,maxerr, CDYN_CONSTRAINT_NORMAL,NULL);
    } 
    else 
    {
        cp=new cDynContactPoint(this,n,r,err,maxerr, CDYN_CONSTRAINT_NORMAL_FRICTION,fr->reference());
        // add constraints tangential to surface
        new cDynContactPoint(this,fr->x(),r,0.0,0.0, CDYN_CONSTRAINT_FRICTION_X,fr->reference());
        new cDynContactPoint(this,fr->y(),r,0.0,0.0, CDYN_CONSTRAINT_FRICTION_Y,fr->reference());
    }

    // add contact point to cache so we can find it if we get it again
    g_dynTheContactPointCache.insert(cp,h);
}

//---------------------------------------------------------------------------

void cDynContact::resetPoints()
{
    assert(type_ == CDYN_CONTACT);
    while (head_ != NULL) 
    {
        delete head_;
    }
    tail_=NULL;
}

//---------------------------------------------------------------------------

void cDynContact::findPoints(const cDynTime& time)
{
    resetPoints();
    if (constraint_ == NULL) 
    {
        prim_->findContacts(time, this);
    } 
    else 
    {
        constraint_->contact(this);
        constraint_->define(time,prim_->nodeA_,prim_->nodeB_);
    }
}

//---------------------------------------------------------------------------
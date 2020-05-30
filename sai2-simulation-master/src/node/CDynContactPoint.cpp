//===========================================================================
/*
    This file is part of the dynamics library.
    Copyright (C) 2014, Artificial Intelligence Laboratory,
    Stanford University. All rights reserved.

    \version   $MAJOR.$MINOR.$RELEASE $Rev: 1242 $
*/
//===========================================================================

//---------------------------------------------------------------------------
#include "node/CDynContact.h"
#include "node/CDynConstraint.h"
//---------------------------------------------------------------------------
cDynContactPoint* cDynContactPoint::free_=NULL;
//---------------------------------------------------------------------------
cDynContactPoint::cDynContactPoint(cDynContact* c, 
    const cDynVector3& n,
    const cDynVector3& r, 
    const double err, 
    const double maxerr, 
    const cDynContactPointType type, 
    cDynFrictionRecord* fr, 
    cDynConstraintInfo* info)
{
    type_=type;
    s_[0]=n;
    r_=r;
    s_[1].crossMultiply(r, n);
    err_=err;
    maxerr_=maxerr;
    fr_=fr;
    info_=info;
    force_=0;
    impulse_ = -100;
    velocity_ = -100;
    deltavelocity_ = -100;
    acceleration_ = -100;

    contact_=c;
    if (c->head_ == NULL) 
    {
        c->head_=c->tail_=this;
        prev_=next_=NULL;
    } 
    else 
    {
        prev_=c->tail_;
        next_=NULL;
        c->tail_->next_=this;
        c->tail_=this;
    }
    c->numPoints_++;
}

//---------------------------------------------------------------------------

cDynContactPoint::~cDynContactPoint()
{
    cDynContact* c=contact_;
    if (c->tail_ == this) c->tail_=prev_;
    if (c->head_ == this) c->head_=next_;
    if (prev_ != NULL) prev_->next_=next_;
    if (next_ != NULL) next_->prev_=prev_;
    c->numPoints_--;
    if (fr_ != NULL) fr_->unreference();
}

//---------------------------------------------------------------------------
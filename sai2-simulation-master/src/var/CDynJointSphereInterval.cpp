//===========================================================================
/*
    This file is part of the dynamics library.
    Copyright (C) 2014, Artificial Intelligence Laboratory,
    Stanford University. All rights reserved.

    \version   $MAJOR.$MINOR.$RELEASE $Rev: 1242 $
*/
//===========================================================================

//---------------------------------------------------------------------------
#include "var/CDynVar.h"
//---------------------------------------------------------------------------
cDynJointSphereInterval* cDynJointSphereInterval::free_=NULL;
//---------------------------------------------------------------------------

cDynJointSphereInterval::cDynJointSphereInterval(cDynJointSphereVar* var,
    const cDynTime t,
    const cDynQuaternion& q, 
    const cDynVector3& v, 
    const cDynVector3& a)
{
    var_=var;
    t_=t;
    q_=q;
    v_=v;
    a_=a;

    if (var->tail_ == NULL) dt_=0.0f;
    else dt_=(double)(t-var->tail_->t_);

    if (var_->tail_ == NULL) 
    {
        var_->head_=var_->tail_=this;
        prev_=next_=NULL;
    } 
    else 
    {
        prev_=var_->tail_;
        next_=NULL;
        var_->tail_->next_=this;
        var_->tail_=this;
    }
}

//---------------------------------------------------------------------------

cDynJointSphereInterval::~cDynJointSphereInterval()
{
    if (var_->tail_ == this) var_->tail_=prev_;
    if (var_->head_ == this) var_->head_=next_;
    if (prev_ != NULL) prev_->next_=next_;
    if (next_ != NULL) next_->prev_=prev_;
}

//---------------------------------------------------------------------------
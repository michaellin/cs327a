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
cDynJointInterval* cDynJointInterval::free_=NULL;
//---------------------------------------------------------------------------

cDynJointInterval::cDynJointInterval(cDynJointVar* var,const cDynTime t,const cDynVector3& q)
{
    assert(var->tail_ == NULL || var->tail_->t_ <= t);
    var_=var;
    t_=t;
    q_=q;
    if (var->tail_ == NULL) dt_=0.0f;
    else dt_=(double)(t-var->tail_->t_);

#ifdef CDYN_EXTENDED
    valid_=0;
#endif

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

cDynJointInterval::~cDynJointInterval()
{
    if (var_->tail_ == this) var_->tail_=prev_;
    if (var_->head_ == this) var_->head_=next_;
    if (prev_ != NULL) prev_->next_=next_;
    if (next_ != NULL) next_->prev_=prev_;
}

//---------------------------------------------------------------------------

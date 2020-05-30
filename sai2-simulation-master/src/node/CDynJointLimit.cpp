//===========================================================================
/*
    This file is part of the dynamics library.
    Copyright (C) 2014, Artificial Intelligence Laboratory,
    Stanford University. All rights reserved.

    \version   $MAJOR.$MINOR.$RELEASE $Rev: 1242 $
*/
//===========================================================================

//---------------------------------------------------------------------------
#include "node/CDynJointLimit.h"
//---------------------------------------------------------------------------

void cDynJointLimitList::insert(cDynJointLimit* limit)
{
    if (limit->list_ != NULL) limit->list_->remove(limit);
    if (limit->list_ == this) return; // already in list
    assert(limit->state_ == -1);
    limit->prev_=tail_;
    limit->next_=NULL;
    if (tail_ != NULL) tail_->next_=limit;
    else head_=limit;
    tail_=limit;
    limit->state_=state_;
    limit->list_=this;
}

//---------------------------------------------------------------------------

void cDynJointLimitList::remove(cDynJointLimit* limit)
{
    assert(limit->list_=this);
    assert(limit->state_ != -1);
    if (limit == tail_) tail_=tail_->prev_;
    if (limit == head_) head_=head_->next_;
    if (limit->next_ != NULL) limit->next_->prev_=limit->prev_;
    if (limit->prev_ != NULL) limit->prev_->next_=limit->next_;
    limit->next_=limit->prev_=NULL;
    limit->state_= -1;
    limit->list_=NULL;
}

//---------------------------------------------------------------------------

void cDynJointLimitList::insert(cDynJoint* joint)
{
    if (joint->bound_[0])
        insert(joint->bound_[0]);
    if (joint->bound_[1])
        insert(joint->bound_[1]);
}

//---------------------------------------------------------------------------

void cDynJointLimitList::pop(cDynJointLimit* limit)
{
    assert(limit->state_ == state_);
    if (limit == head_) return;
    remove(limit);
    // CHECK XXX I think this can just be insert(limit);
    limit->next_=head_;
    limit->prev_=NULL;
    if (head_ != NULL) head_->prev_=limit;
    else tail_=limit;
    head_=limit;
    limit->state_= state_;
    limit->list_=this;
}

//---------------------------------------------------------------------------

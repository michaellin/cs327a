//===========================================================================
/*
    This file is part of the dynamics library.
    Copyright (C) 2014, Artificial Intelligence Laboratory,
    Stanford University. All rights reserved.

    \version   $MAJOR.$MINOR.$RELEASE $Rev: 1242 $
*/
//===========================================================================

//---------------------------------------------------------------------------
#include "node/CDynConstraintList.h"
//---------------------------------------------------------------------------

cDynConstraintList::cDynConstraintList()
{
    head_=tail_=NULL;
}

//---------------------------------------------------------------------------

cDynConstraintList::~cDynConstraintList()
{
    while (head_ != NULL)
        remove(head_);
}

//---------------------------------------------------------------------------

void cDynConstraintList::insert(cDynConstraint* c)
{
    c->next_=NULL;
    c->prev_=tail_;
    if (tail_ != NULL) tail_->next_=c;
    else head_=c;
    tail_=c;
}

//---------------------------------------------------------------------------

void cDynConstraintList::remove(cDynConstraint* c)
{
        if (c == tail_) tail_=tail_->prev_;
        if (c == head_) head_=head_->next_;
        if (c->next_ != NULL) c->next_->prev_=c->prev_;
        if (c->prev_ != NULL) c->prev_->next_=c->next_;
        c->next_=c->prev_=NULL;
}

//---------------------------------------------------------------------------

void cDynConstraintList::check(const cDynTime& lower, const cDynTime& upper)
{
    for (cDynConstraint* c=head_; c != NULL; c=c->next_) 
    {
        c->time_[0]=lower;
        c->time_[1]=upper;
        c->check();
    }
}

//---------------------------------------------------------------------------

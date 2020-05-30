//===========================================================================
/*
    This file is part of the dynamics library.
    Copyright (C) 2014, Artificial Intelligence Laboratory,
    Stanford University. All rights reserved.

    \version   $MAJOR.$MINOR.$RELEASE $Rev: 1242 $
*/
//===========================================================================

//---------------------------------------------------------------------------
#include "node/CDynCollision.h"
//---------------------------------------------------------------------------
//#define HINC(x) ((x+1)%CDYN_BOUND_HISTORY_SIZE)
#define HINC(x) ((x+1) & 0x3)
//---------------------------------------------------------------------------

void cDynBound::update(double value, const cDynTime& tmin, const cDynTime& tmax)
{
    //	value_=value; return;
    int dirty=0;
    
    // clear history that is no longer needed
    while (he_ != hs_ && history_[HINC(he_)].time_ <= tmin) 
    {
        // one of the deleted items is the current so search again
        if  (he_ == index_) dirty++; 
        he_=HINC(he_);
    }

    if (he_ - hs_ == CDYN_BOUND_HISTORY_SIZE) 
    { 
        // list is full so throw away latest entry, but update value
        if (type_ > 0) 
        { 
            // upper
            if (history_[hs_].value_ > value) value=history_[hs_].value_;
        } 
        else 
        {
            if (history_[hs_].value_ < value) value=history_[hs_].value_;
        }
        
        if (hs_ == index_) value_=history_[hs_].value_;
    } 
    else 
    {
        hs_++;
        if (hs_ == CDYN_BOUND_HISTORY_SIZE) hs_=0;
    }
    
    history_[hs_].time_=tmax;
    history_[hs_].value_=value;

    if (type_ > 0) 
    {
        if (history_[hs_].value_ > value_) { value_=history_[hs_].value_; index_=hs_; }
    } 
    else 
    {
        if (history_[hs_].value_ < value_) { value_=history_[hs_].value_; index_=hs_; }
    }

    if (dirty) 
    {
        int c=he_;
        value_=history_[hs_].value_;
        index_=hs_;
        while (c != hs_) 
        {
            if (type_ > 0) 
            {
                if (history_[c].value_ > value_) { value_=history_[c].value_; index_=c; }
            } 
            else 
            {
                if (history_[c].value_ < value_) { value_=history_[c].value_; index_=c; }
            }
            c=HINC(c);
        }
    }
}

//---------------------------------------------------------------------------

void cDynBound::set(double value, cDynTime& time)
{
    he_=hs_=index_=0;
    history_[0].time_=time;
    history_[0].value_=value;
    value_=value;
}

//---------------------------------------------------------------------------

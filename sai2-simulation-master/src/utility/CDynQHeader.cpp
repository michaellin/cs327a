//===========================================================================
/*
    This file is part of the dynamics library.
    Copyright (C) 2014, Artificial Intelligence Laboratory,
    Stanford University. All rights reserved.

    \version   $MAJOR.$MINOR.$RELEASE $Rev: 1242 $
*/
//===========================================================================

//---------------------------------------------------------------------------
#include <cassert>
#include "utility/CDynQHeader.h"
//---------------------------------------------------------------------------
cDynQEvent* cDynQEvent::free_=NULL;
//---------------------------------------------------------------------------

cDynQEvent* cDynQHeader::event(const cDynTime& time,
    void (*f)(const cDynTime& time, void* arg),
    void (*b)(const cDynTime& time, void* arg),
    void* arg)
{
    assert(tail_ == NULL || tail_->time_ <= time);
    cDynQEvent* event=new cDynQEvent(time, f, b, arg);

    // add it to queue
    event->next_=NULL;
    event->prev_=tail_;
    if (event->prev_ == NULL) head_=event;
    else event->prev_->next_=event;
    tail_=event;

    // return pointer
    return(event);
}

//---------------------------------------------------------------------------

cDynQEvent* cDynQHeader::find(const cDynTime& time, void* arg)
{
    cDynQEvent* e=tail_;
    while (e != NULL && e->time_ >= time) 
    {
        if (e->time_ == time && e->arg_ == arg) return(e);
        e=e->prev_;
    }
    return(NULL);
}

//---------------------------------------------------------------------------

void cDynQHeader::remove(cDynQEvent* event)
{
    // remove a event without processing it
    if (event == head_) head_ = event->next_;
    if (event == tail_) tail_ = event->prev_;
    if (event->next_ != NULL) event->next_->prev_=event->prev_;
    if (event->prev_ != NULL) event->prev_->next_=event->next_;
    delete event;
}

//---------------------------------------------------------------------------

void cDynQHeader::advance(const cDynTime& time)
{
    while (head_ != NULL && head_->time_ < time) 
    {
        // disconnect from list;
        // needs to be done first because callback may invoke other calls to advance/backup/remove
        cDynQEvent* event=head_;
        head_=head_->next_;
        if (head_ != NULL) head_->prev_=NULL;
        else tail_ = NULL;

        // process the event
        event->forward();

        // delete the event
        delete event;
    }
}

//---------------------------------------------------------------------------

void cDynQHeader::backup(const cDynTime& time)
{
    while (tail_ != NULL && tail_->time_ >= time) 
    {
        // disconnect from list;
        // needs to be done first because callback may invoke other calls to advance/backup/remove
        cDynQEvent* event=tail_;
        tail_=tail_->prev_;
        if (tail_ != NULL) tail_->next_=NULL;
        else head_=NULL;

        // process the event
        event->backward();

        // delete the event
        delete event;
    }
}

//---------------------------------------------------------------------------
//===========================================================================
/*
    This file is part of the dynamics library.
    Copyright (C) 2014, Artificial Intelligence Laboratory,
    Stanford University. All rights reserved.

    \version   $MAJOR.$MINOR.$RELEASE $Rev: 1242 $
*/
//===========================================================================

//---------------------------------------------------------------------------
#include "utility/CDynLogger.h"
#include "utility/CDynError.h"
#include "utility/CDynTimeEvent.h"
//---------------------------------------------------------------------------

void cDynTimeHeap::Heapify(int i)
{
    int l,r,smallest;

    l=left(i);
    r=right(i);

    if (l <= size_ && less(at(l),at(i)))
    {
        smallest=l;
    }
    else 
    {
        smallest=i;
    }

    if (r <= size_ && less(at(r),at(smallest)))
    {
        smallest=r;
    }
    
    if (smallest != i) 
    {
        // exchange heap[i] and heap[smallest]  
        cDynTimeEvent *tmp;

        tmp=heap_[smallest];
        heap_[smallest]=heap_[i];
        heap_[i]=tmp;

        // update index entries
        heap_[i]->index_=i;
        heap_[smallest]->index_=smallest;
        
        Heapify(smallest);
    }
}

//---------------------------------------------------------------------------

void cDynTimeHeap::insert(cDynTimeEvent* event)
{
    size_++;
    if (size_ == maxsize_) 
    {
        cDynError(CDYN_TIMEHEAP_OVERFLOW);
        return;
    }

    int i=size_;
    while (i > 1 && greater(at(parent(i)),event)) 
    {
        heap_[i]=heap_[parent(i)];
        heap_[i]->index_=i;
        i=parent(i);
    }

    event->heap_=this;
    event->index_=i;
    heap_[i]= event;
}

//---------------------------------------------------------------------------

void cDynTimeHeap::remove(cDynTimeEvent* event)
{
    if (event == active_) active_=NULL;
    int i=event->index_;

    if (i < size_) 
    {
        heap_[i]=heap_[size_];
        heap_[i]->index_=i;
        heap_[size_]=NULL;
        size_--;
        Heapify(i);
    } 
    else 
    {
        heap_[size_]=NULL;
        size_--;
    }

    event->heap_=NULL;
    event->index_=0;
}

//---------------------------------------------------------------------------

void cDynTimeHeap::update(cDynTimeEvent* event)
{
    if (event == NULL) return;
    
    // if updated event is the active event no need to update
    if (event == active_) active_=NULL;

    // move up if needed
    int i=event->index_;
    while (i > 1 && greater(at(parent(i)),event)) 
    {
        heap_[i]=heap_[parent(i)];
        heap_[i]->index_=i;
        i=parent(i);
    }
    
    heap_[i]= event;
    event->index_=i;

    //move down if needed
    Heapify(i);
}

//---------------------------------------------------------------------------

void cDynTimeHeap::next(void)
{
    cDynTimeEvent *event=top();

    if (size_ == 0 || top() == NULL) return;

    // mark active callback event
    active_=event;	

    //make callback
    event->callback();

    // and advance if needed 
    if (active_ == event) 
    {
        switch (event->count_) 
        {
            case 0:  
                event->time_ += event->inc_; 
                break;
            
            case 1:  
                remove(event); 
                return;
            
            default: 
                event->time_ += event->inc_; 
                event->count_--;
                break;
        }
        Heapify(event->index_);
    }
    //debug();
}

//---------------------------------------------------------------------------

void cDynTimeHeap::advance(const cDynTime time)
{
    while (size_ > 0 && top() != NULL && top()->time() <= time) 
    {
        next();
    }
}

//---------------------------------------------------------------------------

void cDynTimeHeap::advance(const int iter)
{
    int count=0;
    while (size_ > 0 && top() != NULL && count < iter) 
    {
        next();
        count++;
    }
}

//---------------------------------------------------------------------------

void cDynTimeHeap::debug(void)
{
    int l,r;

    for (int i=1;i<=size_;i++) 
    {
        if (heap_[i] == NULL) 
        {
            cDynPrintf("position %d has a NULL pointer\n",i);
            return;
        }

        if (heap_[i]->index_ != i) 
        {
            cDynPrintf("position %d index is incorrect %d\n",
                    i,heap_[i]->index_);
        }
        
        r=right(i);
        l=left(i);
        
        if (r < size_ && greater(at(i),at(r))) 
        {
            cDynPrintf("%d right child %d is bigger then parent\n",i,r);
        }
        
        if (l < size_ && greater(at(i),at(l))) 
        {
            cDynPrintf("%d left child %d is bigger then parent\n",i,l);
        }
    }
}

//---------------------------------------------------------------------------
//===========================================================================
/*
    This file is part of the dynamics library.
    Copyright (C) 2014, Artificial Intelligence Laboratory,
    Stanford University. All rights reserved.

    \version   $MAJOR.$MINOR.$RELEASE $Rev: 1242 $
*/
//===========================================================================

//---------------------------------------------------------------------------
#include <math.h>
#include "utility/CDynError.h"
#include "var/CDynState.h"
#include "var/CDynVar.h"
#include "matrix/CDynQuaternion.h"
#include "matrix/CDynVector3.h"
//---------------------------------------------------------------------------
#if defined(WIN32) | defined(WIN64)
#include <float.h>
#define finite _finite
#endif
//---------------------------------------------------------------------------

cDynState::cDynState(const int size): size_(size), n_(0), mark_(0), head_(NULL), tail_(NULL)
{
    x_ = new double[size_];
    dxdt_ = new double[size_];
    mark_=0;

    head_=tail_=NULL;
}

//---------------------------------------------------------------------------

cDynState::~cDynState()
{
    delete[] x_;
    delete[] dxdt_;

    cDynStateEntry *tmp=head_;
    cDynStateEntry *list=head_;

    while (list!=NULL) 
    {
        tmp=list;
        list=list->next_;
        delete tmp;
    }

    head_=tail_=NULL;
}

//---------------------------------------------------------------------------

void cDynState::grow(const int size)
{
    if (size_ >= size) return;

    double* x=new double[size];
    double* dxdt=new double[size];

    for (int i=0;i<n_;i++) 
    {
        x[i]=x_[i];
        dxdt[i]=dxdt_[i];
    }

    delete[] x_;
    delete[] dxdt_;
    x_=x;
    dxdt_=dxdt;
    size_=size;
}

//---------------------------------------------------------------------------

cDynStateEntry* cDynState::insert(cDynVar* var, int s)
{
    if (n_+s >= size_) grow(size_+size_);

    cDynStateEntry* entry=new cDynStateEntry;
    entry->ptr_=var;
    entry->state_=this;
    entry->index_=n_;
    entry->num_=s;
    n_ += s;
    entry->next_=NULL;
    entry->prev_=NULL;
    
    for (int i=0;i<s;i++) 
    {
        x_[entry->index_+i]=0.0f;
        dxdt_[entry->index_+i]=0.0f;
    }

    if (tail_ == NULL) 
    {
        head_=tail_=entry;
    } 
    else 
    {
        entry->prev_=tail_;
        tail_->next_=entry;
        tail_=entry;
    }

    // change of variables mark
    mark(); 

    return(entry);
}

//---------------------------------------------------------------------------

void cDynState::remove(cDynStateEntry* entry)
{
    if (entry == NULL) return;
    if (entry->state_ != this) cDynError(CDYN_ERROR_DESTATE_ERROR1);
    
    //remove from list
    cDynStateEntry* next=entry->next_;
    if (next != NULL) next->prev_=entry->prev_;
    cDynStateEntry* prev=entry->prev_;
    if (prev != NULL) prev->next_=next;
    if (head_ == entry) head_=next;
    if (tail_ == entry) tail_=prev;
    int s=entry->num_;

    //close the gap
    while (next != NULL) 
    {
        next->index_ -= s;
        next=next->next_;
    }

    //update the state vector
    n_ -= s;
    for (int i=entry->index_;i<n_;i++) 
    {
        x_[i]=x_[i+s];
        dxdt_[i]=dxdt_[i+s];
    }

    // change of variables mark
    mark(); 
}

//---------------------------------------------------------------------------

void cDynState::torqueExternalZero()
{
    for (cDynStateEntry* e=head_;e != NULL; e=e->next_) 
    {
        if (e->num_ == 2)
        {
            ((cDynJointVar *)e->ptr_)->torqueExternalZero();
        }
        else
        {
            ((cDynJointSphereVar *)e->ptr_)->torqueExternalZero();
        }
    }
}

//---------------------------------------------------------------------------

void cDynState::update(const bool normalize, const bool err)
{
    // for each variable in state vector update 
    for (cDynStateEntry* e=head_;e != NULL; e=e->next_) 
    {
        // 1 dof joint variable
        if (e->num_ == 2) 
        { 
            if (err) 
            {
                cDynJointVar* var=(cDynJointVar*)e->ptr_;
                var->errorUpdate();
            }
            e->dxdt(0)=e->x(1);  // update dxdt values
        } 

        // spherical joint variable
        else 
        { 
            if (err) 
            {
                cDynJointSphereVar* var=(cDynJointSphereVar*)e->ptr_;
                var->errorUpdate();
            }

            // renormalize quaternion
            if (normalize) 
            {
                double msqr = e->x(0)*e->x(0)
                         + e->x(1)*e->x(1)
                         + e->x(2)*e->x(2)
                         + e->x(3)*e->x(3);
                double mi=1.0f/cDynSqrt(msqr);
                for (int i=0;i<4;i++) e->x(i) *= mi;
            }

            cDynVector3 w;
            w.set(e->x(4),e->x(5),e->x(6));
#if 1
            cDynQuaternion dq;
            dq.velocity(((cDynJointSphereVar *)e->ptr_)->q(), w);
            e->dxdt(0)=dq[0];
            e->dxdt(1)=dq[1];
            e->dxdt(2)=dq[2];
            e->dxdt(3)=dq[3];
#endif
#if 0
            cDynVector3 r;
            r.inversedMultiply(((cDynJointSphereVar *)e->ptr_)->q(), w);
            //r=w;
            
            // note equation also appears in cDynJointSphere, seperated for speed
/*
            e->dxdt(0)=0.5*(-e->x(1)*r[0] - e->x(2)*r[1] - e->x(3)*r[2]);
            e->dxdt(1)=0.5*( e->x(0)*r[0] + e->x(3)*r[1] - e->x(2)*r[2]);
            e->dxdt(2)=0.5*(-e->x(3)*r[0] + e->x(0)*r[1] + e->x(1)*r[2]);
            e->dxdt(3)=0.5*( e->x(2)*r[0] - e->x(1)*r[1] + e->x(0)*r[2]);
*/
            // mirtich
            e->dxdt(0)=0.5*(-e->x(1)*r[0] - e->x(2)*r[1] - e->x(3)*r[2]);
            e->dxdt(1)=0.5*( e->x(0)*r[0] - e->x(3)*r[1] + e->x(2)*r[2]);
            e->dxdt(2)=0.5*( e->x(3)*r[0] + e->x(0)*r[1] - e->x(1)*r[2]);
            e->dxdt(3)=0.5*(-e->x(2)*r[0] + e->x(1)*r[1] + e->x(0)*r[2]);
#endif

        }
    }
}

//---------------------------------------------------------------------------

void cDynState::backup(const cDynTime& time)
{
    // for each variable in state vector backup f-curve
    for (cDynStateEntry* e=head_;e != NULL; e=e->next_) 
    {
        // 1 dof joint variable
        if (e->num_ == 2) 
        { 
            cDynJointVar* var=(cDynJointVar*)e->ptr_;
            var->backup(time);
        } 
        
        // spherical joint variable
        else 
        { 
            cDynJointSphereVar* var=(cDynJointSphereVar*)e->ptr_;
            var->backup(time);
        }
    }
}

//---------------------------------------------------------------------------

void cDynState::update(const cDynTime& time)
{
    // for each variable in state vector update f-curve
    for (cDynStateEntry* e=head_;e != NULL; e=e->next_) 
    {
        // 1 dof joint variable
        if (e->num_ == 2) 
        {
            cDynJointVar* var=(cDynJointVar*)e->ptr_;
            var->update(time);
        } 

        // spherical joint variable
        else 
        { 
            cDynJointSphereVar* var=(cDynJointSphereVar*)e->ptr_;
            var->update(time);
        }
    }
}

//---------------------------------------------------------------------------

void cDynState::advance(const cDynTime& time)
{
    // for each variable in state vector advance the history vector
    for (cDynStateEntry* e=head_;e != NULL; e=e->next_) 
    {
        // 1 dof joint variable
        if (e->num_ == 2) 
        { 
            cDynJointVar* v=(cDynJointVar*)e->ptr_;
            v->advance(time);
        } 
        
        // spherical joint variable
        else 
        { 
            cDynJointSphereVar* v=(cDynJointSphereVar*)e->ptr_;
            v->advance(time);
        }
    }
}

//---------------------------------------------------------------------------

bool cDynState::check()
{
    for (int i=0;i<n_;i++) 
    {
        if (!finite(x_[i]) || x_[i] > 1e6 || x_[i] < -1e6) return(false);
        if (!finite(dxdt_[i]) || dxdt_[i] > 1e6 || dxdt_[i] < -1e6) return(false);
    }
    return(true);
}

//---------------------------------------------------------------------------
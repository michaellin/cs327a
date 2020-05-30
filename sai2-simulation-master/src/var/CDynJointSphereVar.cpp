//===========================================================================
/*
    This file is part of the dynamics library.
    Copyright (C) 2014, Artificial Intelligence Laboratory,
    Stanford University. All rights reserved.

    \version   $MAJOR.$MINOR.$RELEASE $Rev: 1242 $
*/
//===========================================================================

//---------------------------------------------------------------------------
#include "utility/CDynError.h"
#include "utility/CDynLogger.h"
#include "var/CDynInterpolate.h"
#include "var/CDynVar.h"
#include "matrix/CDynFrame.h"
#include "matrix/CDynQuaternion.h"
//---------------------------------------------------------------------------
cDynFrame DeDebugFrame;
//---------------------------------------------------------------------------

void cDynJointSphereVar::state(cDynState* s)
{
    cDynQuaternion qh; qh=qCurrent();
    cDynVector3 vh; vh=vCurrent();
    cDynVector3 ah; ah=aCurrent();
    if (state_ != NULL) 
    { 
        // if already part of some other state delete it
        //save last position, just in case
        qh=q();	
        vh=v();
        ah=a();
        delete state_;
        state_=NULL;
    }
    if (s) 
    {
        state_=s->insert(this);
        q(qh);
        v(vh);
        a(ah);
    } 
    else 
    {
        qCurrent(qh);
        vCurrent(vh);
        aCurrent(ah);
    }
}

//---------------------------------------------------------------------------

void cDynJointSphereVar::backup(const cDynTime t)
{
    if (tail_ != NULL && tail_->t_ > t) 
    { 
        // need to go backward
        while (tail_->prev_ != NULL && tail_->prev_->t_ > t)
        {
            delete tail_;
        }

        if (tail_->prev_ == NULL) 
        {
            // This case should not happen since we are
            // backing up past the end of time, but just
            // in case do the right thing.
            tail_->t_=t;
        } 
        else if (tail_->prev_->t_ == t) 
        {
            delete tail_;
        } 
        else 
        {
            double pt= (double)(t - tail_->prev_->t_);

            double u=pt/tail_->dt_;
            for (int i=0;i<4;i++)
            {
                tail_->q_[i]= tail_->prev_->q_[i] + u*(tail_->q_[i] - tail_->prev_->q_[i]);
            }
            tail_->q_.normalize();
            for (int i=0;i<3;i++)
            {
                tail_->v_[i]= tail_->prev_->v_[i] + u*(tail_->v_[i] - tail_->prev_->v_[i]);
            }
            for (int i=0;i<3;i++)
            {
                tail_->a_[i]= tail_->prev_->a_[i] + u*(tail_->a_[i] - tail_->prev_->a_[i]);
            }
            tail_->t_=t;
            tail_->dt_=pt;
        }
    }
    q_=tail_->q_;
    v_=tail_->v_;
    a_=tail_->a_;
    t_=t;
    q(q_);
    v(v_);
    a(a_);
}

//---------------------------------------------------------------------------

void cDynJointSphereVar::errorUpdate()
{
    if (errorExists_) 
    {
        cDynVector3 s;
        cDynQuaternion err;
        const cDynQuaternion& r=q();
        s.multiply(r,error_);
        err.velocity(r,s);
        err += r;
        q(err);
/*
        cDynPrintf("err: [%5.5e, %5.5e, %5.5e ]\n",
            error_[0],
            error_[1],
            error_[2]
        );
        cDynQuaternion x2(cDynVector3(1.0f,0.0f,0.0f),10.0*error_[0]);
        cDynQuaternion y2(cDynVector3(0.0f,1.0f,0.0f),10.0*error_[1]);
        cDynQuaternion z2(cDynVector3(0.0f,0.0f,1.0f),10.0*error_[2]);
        DeDebugFrame = x2*y2*z2;
*/
        error_.zero();
        errorExists_=false;
    }
}

//---------------------------------------------------------------------------

void cDynJointSphereVar::update(const cDynTime t)
{
    // add new values
    q_=q();
    v_=v();
    a_=a();
    
    new cDynJointSphereInterval(this,t,q_,v_,a_);
    t_=t;
}

//---------------------------------------------------------------------------

void cDynJointSphereVar::advance(const cDynTime t)
{
    if (head_==NULL) return;

    cDynJointSphereInterval* cutoff=head_;
    while (cutoff->next_ != NULL && cutoff->next_->t_ <= t)
        cutoff=cutoff->next_;
    while (head_ != cutoff)
        delete head_;
}

//---------------------------------------------------------------------------

void cDynJointSphereVar::lookup(const cDynTime t) const
{
    cDynJointSphereInterval* p=tail_;
    cDynJointSphereInterval* o=NULL;
    
    while (p != NULL && p->t_ > t) 
    {
        o=p;
        p=p->prev_;
    }
    
    if (p == NULL) 
    {
        if (o == NULL)
        {
            // lookup below min time
            cDynError(CDYN_ERROR_FIX);	
        }
        else 
        {  
            // done to avoid problems with round off errors
            q_=o->q_;
            v_=o->v_;
            a_=o->a_;
            t_= t;
            return;
        }
    }
    
    if (o == NULL) 
    { 
        // lookup beyond max time
        // done to avoid problems with round off
        q_=p->q_;
        v_=p->v_;
        a_=p->a_;
        t_=t; return;
    }

    if (t == o->t_) 
    {
        q_=o->q_; v_=o->v_; a_=o->a_;
        t_= t; return;
    }
    
    if (t == p->t_) 
    {
        q_=p->q_; v_=p->v_; a_=p->a_;
        t_= t;
        return;
    }

    // otherwise linearly interpolate the result
    double pt= (double)(t - p->t_);

    double u=pt/o->dt_;
    for (int i=0;i<4;i++)
    {
        q_[i]= p->q_[i] + u*(o->q_[i] - p->q_[i]);
    }
    q_.normalize();
    
    for (int i=0;i<3;i++)
    {
        v_[i]= p->v_[i] + u*(o->v_[i] - p->v_[i]);
    }
    for (int i=0;i<3;i++)
    {
        a_[i]= p->a_[i] + u*(o->a_[i] - p->a_[i]);
    }
    t_=t;
}

//---------------------------------------------------------------------------

size_t cDynJointSphereVar::save(void* ptr, size_t size)
{
    char* p=(char *)ptr;
    int n=0;
    for (cDynJointSphereInterval* q=head_;q != NULL; q=q->next_) n++;

    size_t s=sizeof(int) + n*(sizeof(cDynTime) + (4+3+3)*sizeof(double));

    if (size >= s) 
    {
        *((int *)p) = n; p += sizeof(int);
        for (cDynJointSphereInterval* q=head_;q != NULL; q=q->next_) 
        {
            *((cDynTime *)p) = q->t_; p += sizeof(cDynTime);
            {for (int i=0;i<4;i++) 
            {
                *((cDynTime *)p) = q->q_[i]; p += sizeof(double);
            }}
            {for (int i=0;i<3;i++) 
            {
                *((cDynTime *)p) = q->v_[i]; p += sizeof(double);
            }}
            {for (int i=0;i<3;i++) 
            {
                *((cDynTime *)p) = q->a_[i]; p += sizeof(double);
            }}
        }
    }

    return(s);
}

//---------------------------------------------------------------------------

size_t cDynJointSphereVar::restore(void* ptr, size_t size, cDynTime offset)
{
    char* p=(char *)ptr;

    if (size < sizeof(int)) return(0);
    int n= *((int *)p); p += sizeof(int);

    size_t s=sizeof(int) + n*(sizeof(cDynTime) + (4+3+3)*sizeof(double));
    if (s > size) return(0);

    // delete all information currently being stored
    while (head_ != NULL)
        delete head_;

    for (int j=0;j<n;j++) 
    {
        t_ = *((cDynTime *)p);  p += sizeof(cDynTime);
        {for (int i=0;i<4;i++) 
        {
            q_[i] = *((cDynTime *)p); p += sizeof(double);
        }}
        {for (int i=0;i<3;i++) 
        {
            v_[i] = *((cDynTime *)p); p += sizeof(double);
        }}
        {for (int i=0;i<3;i++) 
        {
            a_[i] = *((cDynTime *)p); p += sizeof(double);
        }}
        t_ += offset;
        new cDynJointSphereInterval(this,t_,q_,v_,a_);
    }
    
    q(q_); v(v_); a(a_);

    return(s);
}

//---------------------------------------------------------------------------

#ifdef CDYN_DEBUG
void cDynJointSphereVar::display(const char* str)
{
    cDynPrintf ("%s = \n",str);
    for (cDynJointSphereInterval* q=head_;q != NULL; q=q->next_) 
    {
        cDynPrintf("%5.9f %5.9f %5.9f %5.9f %5.9f\n",
            q->t_, q->q_[0], q->q_[1], q->q_[2], q->q_[3]
        );
    }
}

//---------------------------------------------------------------------------
#endif
//---------------------------------------------------------------------------
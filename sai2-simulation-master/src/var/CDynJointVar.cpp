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
//---------------------------------------------------------------------------
#ifdef CDYN_EXTENDED
static cDynInterpolate interpolate;
#endif
//---------------------------------------------------------------------------

void cDynJointVar::state(cDynState* s)
{
    double qh=qCurrent();
    double vh=vCurrent();
    double ah=aCurrent();
    if (state_ != NULL) 
    {	
        // if already part of some other state delete it
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

void cDynJointVar::backup(const cDynTime t)
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
#ifdef CDYN_EXTENDED
            if (!tail_->valid_) 
            {
                interpolate.generate(tail_->dt_, tail_->a_,tail_->prev_->q_,tail_->q_);
                tail_->valid_++;
            }
            tail_->t_=t;
            interpolate.solve(pt,tail_->a_,tail_->q_);
            tail_->dt_=pt;
            tail_->valid_=0;
#else
            double u=pt/tail_->dt_;
            for (int i=0;i<3;i++)
            {
                tail_->q_[i]= tail_->prev_->q_[i] + u*(tail_->q_[i] - tail_->prev_->q_[i]);
            }

            tail_->t_=t;
            tail_->dt_=pt;
#endif
        }
    }
    // update state with previous values
    q(tail_->q_[0]); v(tail_->q_[1]); a(tail_->q_[2]);
    t_=t;
    q_[0]=q();
    q_[1]=v();
    q_[2]=a();
}

//---------------------------------------------------------------------------

void cDynJointVar::update(const cDynTime t)
{
    t_=t;
    q_[0]=q();
    q_[1]=v();
    q_[2]=a();
    new cDynJointInterval(this,t,q_);
}

//---------------------------------------------------------------------------

void cDynJointVar::advance(const cDynTime t)
{
    if (head_==NULL) return;

    cDynJointInterval* cutoff=head_;
    
    while (cutoff->next_ != NULL && cutoff->next_->t_ <= t)
        cutoff=cutoff->next_;

    while (head_ != cutoff)
        delete head_;
}

//---------------------------------------------------------------------------

void cDynJointVar::lookup(const cDynTime t)
{
    cDynJointInterval* p=tail_;
    cDynJointInterval* o=NULL;

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
            // done to avoid problems with round off
            q_=o->q_; t_= t; return;
        }
    }

    if (o == NULL) 
    { 
        // lookup beyond max time
        // done to avoid problems with round off
        q_=p->q_; t_= t; return;
    }

    // check end cases to see if we can save a little time
    if (t == o->t_) { q_=o->q_; t_= t; return; }
    if (t == p->t_) { q_=p->q_; t_= t; return; }

    double pt= (double)(t - p->t_);
#ifdef CDYN_EXTENDED
    if (!p->valid_) {
        interpolate.generate(p->dt_, p->a_,p->prev_->q_,p->q_);
        p->valid_++;
    }
    t_=t;
    interpolate.solve(pt,p->a_,q_);
#else
    double u=pt/o->dt_;
    for (int i=0;i<3;i++)
    {
        q_[i]= p->q_[i] + u*(o->q_[i] - p->q_[i]);
    }
    t_=t;
#endif
}

//---------------------------------------------------------------------------

cDynTime cDynJointVar::lastEqual(const double value, const double err, const cDynTime time, const int dir) const
{
    // find last (latest) point in time before time=time
    // where joint variable is = value +- err.
    // it is required that if:
    //    dir == -1 : position at time=time < value - err
    //    dir ==  1 : position at time=time > value + err

    cDynJointInterval* p=tail_;
    if (p == NULL) 
    { 
        // something bogus empty f-curve
        return(CDYN_TIME_UNINITIALIZED);
    }
    
    cDynJointInterval* o=p;
    p=p->prev_;

    while (p != NULL && p->t_ > time) 
    {
        o=p;
        p=p->prev_;
    }
    
    double bound=(dir == -1)?value-err:value+err;
    if (dir == -1) 
    {
        while (p != NULL && p->q_[0] < bound) 
        {
            o=p; p=p->prev_;
        }
    } 
    else 
    {
        while (p != NULL && p->q_[0] > bound) 
        {
            o=p; p=p->prev_;
        }
    }

    if (p == NULL) 
    {
        if (o->t_ <= time && fabs(o->q_[0] - value) <= err) return(o->t_);
        else return(CDYN_TIME_UNINITIALIZED);
    }
/*
    if (dir == -1) 
        assert(o->q_[0] <= bound);
    else
        assert(o->q_[0] >= bound);
*/
    double dop=o->q_[0] - p->q_[0];
    if (dop < 1e-8 && dop > -1e-8) 
    {
        if (time < o->t_) return(time);
        else return(o->t_); // curve is flat, return lastest point
    }
    
    cDynTime dt=o->t_ - p->t_;
    double dvp=bound - p->q_[0];
    cDynTime t=(dvp*dt)/dop;
    assert (t > -1e-8 && t < dt + 1e-8);
    return(p->t_+t);
}

//---------------------------------------------------------------------------

size_t cDynJointVar::save(void* ptr, size_t size)
{
    char* p=(char *)ptr;
    int n=0;
    for (cDynJointInterval* q=head_;q != NULL; q=q->next_) n++;

    size_t s=sizeof(int) + n*(sizeof(cDynTime) + 3*sizeof(double));

    if (size >= s) 
    {
        *((int *)p) = n; p += sizeof(int);
        for (cDynJointInterval* q=head_;q != NULL; q=q->next_) 
        {
            *((cDynTime *)p) = q->t_; p += sizeof(cDynTime);
            for (int i=0;i<3;i++) 
            {
                *((cDynTime *)p) = q->q_[i]; p += sizeof(double);
            }
        }
    }
    return(s);
}

//---------------------------------------------------------------------------

size_t cDynJointVar::restore(void* ptr, size_t size, cDynTime offset)
{
    char* p=(char *)ptr;

    if (size < sizeof(int)) return(0);
    int n= *((int *)p); p += sizeof(int);

    size_t s=sizeof(int) + n*(sizeof(cDynTime) + 3*sizeof(double));
    if (s > size) return(0);

    // delete all information currently being stored
    while (head_ != NULL)
        delete head_;

    for (int j=0;j<n;j++) 
    {
        t_ = *((cDynTime *)p);  p += sizeof(cDynTime);
        for (int i=0;i<3;i++) 
        {
            q_[i] = *((cDynTime *)p); p += sizeof(double);
        }
        t_ += offset;
        new cDynJointInterval(this,t_,q_);
    }
    q(q_[0]); v(q_[1]); a(q_[2]);
    return(s);
}

//---------------------------------------------------------------------------

#ifdef CDYN_DEBUG
void cDynJointVar::display(const char* str)
{
    cDynPrintf ("%s = \n",str);
    for (cDynJointInterval* q=head_;q != NULL; q=q->next_) {
        cDynPrintf("%5.9f %5.9f\n",
            q->t_, q->q_[0]
        );
    }
}

//---------------------------------------------------------------------------
#endif
//---------------------------------------------------------------------------
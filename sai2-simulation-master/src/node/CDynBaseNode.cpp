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
#include "matrix/CDynQuaternion.h"
#include "node/CDynBaseNode.h"
#include "node/CDynJointLimit.h"
#include "object/CDynObject.h"
#include "node/CDynWorld.h"
#include "dynamics/CDynDynamics.h"
#include "utility/CDynError.h"
//---------------------------------------------------------------------------

void cDynBaseNode::link(cDynObject* obj)
{
    // assert the obj is not connected somewhere else
    assert(obj != NULL && obj->baseNode() == NULL);

    // assert that baseNode not already used by other tree
    assert(root_ == NULL);
    assert(world() != NULL);
    
    beginChange();
        // homeFrame for root is the current world frame
        obj->homeFrame()=world()->frame.top();
        obj->globalFrame()=world()->frame.top();
        root_=obj;
        // rootNode_=((cDynNode *)obj)->root();
    endChange();
}

//---------------------------------------------------------------------------

void    cDynBaseNode::unlink()
{
    collision_=NULL;
    root_->baseSet(NULL,true);
}

//---------------------------------------------------------------------------

cDynBaseNode::cDynBaseNode(cDynWorld* world)
{
    world_=world;
    model_=CDYN_FULL;
    status_=CDYN_INACTIVE;

    // allocate state space for the base
    state_=new cDynState();

    // allocate the contact list for the base
    contact_= new cDynContactList(this);

    collision_=NULL;

    // allocate space for joint limit list
    limit_[CDYN_LIMIT_ACTIVE]=new cDynJointLimitList(CDYN_LIMIT_ACTIVE);
    limit_[CDYN_LIMIT_INVALID]=new cDynJointLimitList(CDYN_LIMIT_INVALID);

    // create tree of nodes
    // rootNode_=NULL;
    root_=NULL;
    
    //get default gravity
    gravity=world_->gravity;

    prev_=next_=NULL;

    fwdInCallback_=NULL;
    fwdInSphereCallback_=NULL;
    fwdOutCallback_=NULL;
    fwdOutSphereCallback_=NULL;

    // create a integration event for base
    updateEvent_=new cDynTimeEvent(world->time(), 1, cDynTime(1.0/120.0), CDYN_TIMEEVENT_FOREVER, &cDynIntegrateCallback, (void *)this);

    // set integrator: NOTE TO BE MOVED TO ITS OWN CLASS
    integrator_=CDYN_EULER_EXPLICIT;
    iterNumMax_=1;
    iterErrorMax_=0.05f;  // relative error bound
    discontinuity_=true;
    
    callnum_=1;
    queryTime_ = -1.0f;
    dynamicsTime_ = -1.0f;
}

//---------------------------------------------------------------------------

void cDynBaseNode::insert(cDynBaseNode** head, cDynBaseNode** tail)
{
    if (*tail == NULL) 
    {
        prev_=next_=NULL;
        *head=*tail=this;
    } 
    else 
    {
        prev_=*tail;
        next_=NULL;
        (*tail)->next_=this;
        *tail=this;
    }
}

//---------------------------------------------------------------------------

void cDynBaseNode::remove(cDynBaseNode** head, cDynBaseNode** tail)
{
    //disconnect from parents children list
    if (world_) 
    {
        if (*tail == this) *tail=prev_;
        if (*head == this) *head=next_;
        if (prev_ != NULL) prev_->next_=next_;
        if (next_ != NULL) next_->prev_=prev_;
        world_=NULL;
    }
}

//---------------------------------------------------------------------------

cDynBaseNode::~cDynBaseNode()
{
    //remove all joint variables
    //reset();

    status(CDYN_INACTIVE);

    if (root_)
        unlink();

    if (world_)
        world_->remove(this);
}

//---------------------------------------------------------------------------

void cDynBaseNode::beginChange()
{
    if (world()) 
    {
        backup(world()->time());
    }
}

//---------------------------------------------------------------------------

void cDynBaseNode::endChange()
{
    root()->baseSet(this,true);
        cDynamicsInitialize(root(), callnum_);
    if (world()) 
    {
        world()->collision().checkObjects(world_->time(),CDYN_TESTOPTION_IGNORE,this);
        world()->collision().reset();
    }
}

//---------------------------------------------------------------------------

void cDynBaseNode::status(cDynStatusType status)
{
    if (status == status_) return;
    switch (status) 
    {
        case CDYN_ISOLATED:
            // add integrate event to world if needed
            if (!updateEvent_->attached()) 
            {
                updateEvent_->time(world_->time());
                world_->insert(updateEvent_);
            }
            beginChange();
                collision_= NULL;
                contact()->updatePersistent(world_->time());
            endChange();
            break;

        case CDYN_ACTIVE:
            // add integrate event to world if needed
            if (!updateEvent_->attached()) 
            {
                updateEvent_->time(world_->time());
                world_->insert(updateEvent_);
            }
            beginChange();
                collision_= &(world_->collision());
                contact()->updatePersistent(world_->time());
            endChange();
            break;

        case CDYN_INACTIVE:
            // remove integrate event from world
            if (updateEvent_->attached()) 
            {
                world_->remove(updateEvent_);
            }
            beginChange();
                collision_= NULL;
                contact()->updatePersistent(world_->time());
            endChange();
            break;

        case CDYN_SLEEP:
            // remove integrate event from world
            if (updateEvent_->attached()) 
            {
                world_->remove(updateEvent_);
            }
            beginChange();
                collision_= &(world_->collision());
                contact()->updatePersistent(world_->time());
            endChange();
            break;
    }

    status_=status;
    dynamicsTime_=queryTime_=world_->time();
    callnum_=0;
    cDynamicsFwdDynamics(root(), &gravity, callnum_);
    state_->update(world_->time());
}

//---------------------------------------------------------------------------

double cDynBaseNode::potentialEnergy(const cDynTime& time)
{
    fwdDynamics(time, false);
    double PE=cDynamicsPotentialEnergy(root(),&gravity);
    return(0.5*PE);
}

//---------------------------------------------------------------------------

double cDynBaseNode::kineticEnergy(const cDynTime& time)
{
    fwdDynamics(time, false);
    double KE=cDynamicsKineticEnergy(root());
    return(0.5*KE);
}

//---------------------------------------------------------------------------

void cDynBaseNode::backup(const cDynTime time)
{
    if (time < updateEvent_->time()) 
    {
        // cDynPrintf("backup %s %13.9f %13.9f  %13.9f\n", data(), time, updateEvent_->time(),dynamicsTime_);
        rollback_.backup(time);
        if (updateEvent_->time() > time) 
        {
            // change state variables to older state
            state_->backup(time);      

            // set timer to integrate time from new state
            updateEvent_->time(time);  
        }
        // flush contact list if later than backup time
        // XXX: fucked up!
        if (contact()->time() < time)
        {
            // NEED TO CHECK
            //contact()->flush();
            contact()->updatePersistent(time);
        }

        dynamicsTime_= -1.0;
        callinc();
    }
}

//---------------------------------------------------------------------------

void cDynBaseNode::fwdDynamics(const cDynTime& time,bool config)
{
    if (dynamicsTime_ != time) 
    {
        callnum_=0;
        queryTime_=time;
        dynamicsTime_=time;
        if (config)
            cDynamicsInitFwdDynamicsConfig(root(), callnum_);
        else
            cDynamicsFwdDynamics(root(), &gravity, callnum_);
    }
}

//---------------------------------------------------------------------------

void cDynBaseNode::integrate(const cDynTime& time,cDynTime& inc)
{
    static int flag = 0;
    double* x=state_->x_;
    double* dxdt=state_->dxdt_; 
    int n=state_->n_;
    double it=(double)inc;
    cDynTime  ntime=time+inc;

    if (n == 0) 
    {
        updateEvent_->time(ntime);
        return;
    }

    // recompute the dynamics if some information is dirty
    if (dynamicsTime_ != time || discontinuity_) 
    {
        //printf("integate BGN %5.5f %s\n",time,name()->string());
        callnum_=0;
        queryTime_=time;
        dynamicsTime_=time;
        cDynamicsFwdDynamics(root(), &gravity, callnum_);
    }

    if (discontinuity_) 
    {
        state_->update(false,false);
        state_->update(time);
        discontinuity_=false;
    }

    bool err =contact()->error();

    switch(integrator_) 
    {
        //-----------------------------------------------------------------------
        case CDYN_EULER_EXPLICIT: 
        //-----------------------------------------------------------------------
        { 
            //
            // forward
            // Xn+1 = Xn + h*d(Xn)dt
            //
            for (int i=0;i<n;i++) 
            {
                x[i] += it*dxdt[i];
            }
        }
        break;

        //-----------------------------------------------------------------------
        case CDYN_EULER_MODIFIED: 
        //-----------------------------------------------------------------------
        { 
            //
            // midpoint
            // Xn+1 = Xn + h*d(Xn+0.5)dt where Xn+0.5 = Xn + 0.5*h*d(Xn)dt
            //

            double *xs=(double *)CDYN_ALLOCA(n*sizeof(double));

            {for (int i=0;i<n;i++)
                xs[i]=x[i];
            }

            double it_2 = 0.5f*it;
            {for (int i=0;i<n;i++) 
                x[i] += it_2*dxdt[i];
            }

            state_->update(true,false);
            cDynamicsFwdDynamics(root(), &gravity, callnum_);
            {for (int i=0;i<n;i++)
                x[i] = xs[i] + it*dxdt[i];
            }
        }
        break;

        //-----------------------------------------------------------------------
        case CDYN_EULER_HEUN: 
        //-----------------------------------------------------------------------
        { 
            //
            // trapezoidal: iterative Heun Method:p596 in Numerical Method for Engineers 2nd by Chapra
            // Xn+1 = Xn + 0.5*h*(d(Xn)dt + d(Xn+1)dt)
            //

            double *xs=(double *)CDYN_ALLOCA(n*sizeof(double));
            double *xe=(double *)CDYN_ALLOCA(n*sizeof(double));
            double *dxdts=(double *)CDYN_ALLOCA(n*sizeof(double));

            {for (int i=0;i<n;i++)
            {
                xs[i]=x[i];
                dxdts[i]=dxdt[i];
                x[i] += it*dxdt[i];	  // initial guess : use Explicit Euler

            }}

            int errorIndex = 0;
            int iter=0;

            while (iter < iterNumMax_ && errorIndex >= 0)
            {
                state_->update(true,false);
                cDynamicsFwdDynamics(root(), &gravity, callnum_);
                {for (int i=0;i<n;i++) 
                {
                    dxdt[i] = 0.5f*(dxdts[i] + dxdt[i]);
                    xe[i] = xs[i] + it*dxdt[i];
                }}
#ifndef CDYN_DEBUG
                if (iter < iterNumMax_ - 1)  // don't check error for the last iteration
#endif 
                    errorIndex = ErrorRelative(xe,x,n);
#ifdef CDYN_DEBUG
                if (iter == iterNumMax_-1 && errorIndex>=0)
                    cDynPrintf("++++++++++ rel error at %d = %f\n",errorIndex,fabs((xe[errorIndex]-x[errorIndex])/xe[errorIndex]));
#endif 
                {for (int i=0;i<n;i++)
                    x[i]=xe[i];
                }
                iter++;
            }
#ifdef CDYN_DEBUG
            if (iterNumMax_ > 1)
            {
                if (iter>2)
                    cDynPrintf("########## # of iteration = %d\n",iter);
                if (errorIndex>=0)
                    cDynPrintf(">>>>>>>>>> nonconvergence at i=%d\n",errorIndex);
            }
#endif
        }
        break;

        //-----------------------------------------------------------------------
        case CDYN_EULER_IMPLICIT: 
        //-----------------------------------------------------------------------    
        { 
            //
            // backward : p152,p158 in Numerical Method for Engineers 2nd by Chapra
            // p58,p71,p413 in An Intro. to Numerical Analysis by Atkinson
            // Xn+1 = Xn + h*d(Xn+1)dt
            //

            double *xo=(double *)CDYN_ALLOCA(n*sizeof(double));
            double *xs=(double *)CDYN_ALLOCA(n*sizeof(double));
            double *xe=(double *)CDYN_ALLOCA(n*sizeof(double));
            double *gs=(double *)CDYN_ALLOCA(n*sizeof(double));
            double *g=(double *)CDYN_ALLOCA(n*sizeof(double));

            {for (int i=0;i<n;i++)
            {  
                xs[i]=xo[i]=x[i];	    // save Xn == xo  // stays constant
                                        // save Xn^i-1 == xs and g(Xn^i-1) == fs
                gs[i] = -it*dxdt[i];    // g(xs) = Xs - Xo - it*d(Xs)dt = -it*dxdt;
                x[i] -= gs[i];          // Xn^i == x: initial guess : use Explicit Euler
            }}

            int errorIndex = 0;
            int iter=0;

            while (iter < iterNumMax_ && errorIndex >= 0)
            {
                state_->update(true,false);
                cDynamicsFwdDynamics(root(), &gravity, callnum_);  // compute d(Xn^i)dt == dxdt
                {for (int i=0;i<n;i++) 
                {
                    //
                    // g(Xn) = Xn - Xo - it*d(Xn)dt = 0 // find a root Xn, note that Xo stays constant.
                    // Xn^i+1 = Xn^i - g(Xn^i)/dg(Xn^i)dx : Newton-Raphson method
                    // Secant method: approx of Newton-Raphson method
                    // df(Xn^i)dx = (g(Xn^i-1) - g(Xn^i))/(Xn^i-1 - Xn^i)  : Secant method : approximation
                    //

                    g[i] = x[i] - xo[i] - it*dxdt[i];
                    if (fabs(g[i])<=iterErrorMax2_ 
                        || fabs(g[i]-gs[i])<=iterErrorMax2_)
                    //	|| fabs(x[i]-xs[i])<=iterErrorMax2_)
                        xe[i] = x[i];
                    else
                        xe[i] = x[i] - g[i]*(x[i]-xs[i])/(g[i]-gs[i]);  // compute Xn^i+1 == xe
                }
                }
#ifndef CDYN_DEBUG
                if (iter < iterNumMax_ - 1)  // don't check error for the last iteration
#endif 
                    errorIndex = ErrorRelative(xe,x,n);
#ifdef CDYN_DEBUG
                if (iter == iterNumMax_-1 && errorIndex>=0)
                    cDynPrintf("++++++++++ rel error at %d = %f\n",errorIndex,fabs((xe[errorIndex]-x[errorIndex])/xe[errorIndex]));
#endif 
                {for (int i=0;i<n;i++) 
                {
                    gs[i]=g[i];			// save g(Xn^i-1) == gs
                    xs[i]=x[i];			// save Xn^i-1 == xs 
                    x[i]=xe[i];			// save Xn^i = x
                }}
                iter++;
            }
#ifdef CDYN_DEBUG
            if (iterNumMax_ > 1)
            {
                if (iter>2)
                    cDynPrintf("########## # of iteration = %d\n",iter);
                if (errorIndex>=0)
                    cDynPrintf(">>>>>>>>>> nonconvergence at i=%d\n",errorIndex);
            }
#endif
        }
        break;
    }
    state_->update(true,err);
    
    if (integrator_ == CDYN_EULER_EXPLICIT) 
    {
        dynamicsTime_=time;
        queryTime_=time;
        callnum_=0;
    } 
    else 
    {
        dynamicsTime_= -1;
        queryTime_= -1;
        callnum_=0;
    }
/*
    dynamicsTime_=ntime;
    queryTime_=ntime;
    callnum_=0;
    //printf("integate END %5.5f %s\n",ntime,name()->string());
    cDynamicsFwdDynamics(rootNode_,&gravity);
*/
    contact()->error(false);

    state_->update(ntime);
    updateEvent_->time(ntime);
    //world_->collisionSet(ntime);

    if (false&&flag) 
    {
        //world_->draw(ntime);
        cDynPrintf("integrate end %5.9f\n",time+inc);

        double PE=cDynamicsPotentialEnergy(root(),&gravity);
        double KE=cDynamicsKineticEnergy(root());

        cDynPrintf("--> ENERGY %s %5.9f KE=%5.9f PE=%5.9f\n",
            data(),
            KE+PE, KE, PE
        );

    }
}

//---------------------------------------------------------------------------
/*!
    absError = xe - x;
    relError = abs(absError/xe) 
    if (xe != 0) relError < eps
    else absError < errorBound(???)
    errorBound = 0.5 * 10^(-m) where m is the number of significant digits 
    of x with respect to xe
    if |xe-x|<=e*|xe|, okay!
    OPTIMIZE THIS such as xe(1-e)<=x okay!
*/
//---------------------------------------------------------------------------

int cDynBaseNode::ErrorRelative(double *xe,double *x,int n)
{
    for (int i=0;i<n;i++)
    {
        //
        // if |xe-x|<=e*e, okay!
        // if |xe-x|>e*e
        //		if |xe|<=e, return i  since |xe-x|>e*e>=e*|xe|
        //		else if |xe-x|>e*|xe|, return i
        //
        double absError = (xe[i]>x[i]) ? xe[i]-x[i] : x[i]-xe[i];
        
        // errorBound == eps*eps : working okay!
        if (absError > iterErrorMax2_) 
        {
            double absXe = (xe[i]>0.0f) ? xe[i] : -xe[i];
            if (absXe>iterErrorMax_)
            {
                //  |xe-x| > e*|xe| with |xe-x|>e*e so, if e>=|xe| return i
                //										else if |xe-x|>e*|xe|, return i
                if (absError > absXe*iterErrorMax_)
                    return i;	
            }
            else
                return i; // |xe|<=e, return i  since |xe-x|>e*e>=e*|xe|
        }
    }
    return -1;
}

//---------------------------------------------------------------------------

// this is the callback function for a backup event
static void cDynBaseNodeBackup(const cDynTime& time, void* arg)
{
    cDynBaseNode* base=(cDynBaseNode *)arg;

    base->backup(time);
}

//---------------------------------------------------------------------------

void cDynBaseNode::dependent(cDynBaseNode* base)
{
    // add events in queue to force rollback of base (base) if (this) base
    // is backed up
    cDynTime t=contact()->time();
    assert( t == base->contact()->time());
    if (!rollback_.find(t,(void *)base)) 
    {
        rollback_.event(t,NULL,&cDynBaseNodeBackup,(void *)base);
        base->rollback_.event(t,NULL,&cDynBaseNodeBackup,(void *)this);
    }
}

//---------------------------------------------------------------------------

void cDynBaseNode::checkInvalids(cDynTime& t)
{
    cDynJointLimit* l=limit_[CDYN_LIMIT_INVALID]->next();
    while (l != NULL) 
    {
        if (l->check(t) >= 0) 
        {
            cDynJointLimit* tmp=l;
            l=limit_[CDYN_LIMIT_INVALID]->next(l);

            limit_[CDYN_LIMIT_INVALID]->remove(tmp);
            limit_[CDYN_LIMIT_ACTIVE]->insert(tmp);

#if 0
            cDynPrintf("Joint %s marked as valid\n",tmp->joint()->data());
#endif
        } 
        else 
            l=limit_[CDYN_LIMIT_INVALID]->next(l);
    }
}
#if 0
const cDynTime cDynBaseNode::checkJointLimits(const cDynTime& time)
{
    if (status_ != CDYN_ACTIVE) return(time);
    static int flag=0;
    flag=(cDynDebugLevel > 3)?1:0;
    cDynTime t=time;
    cDynTime tmax=time;
    cDynTime tmin=world()->time();

    cDynJointLimitList* list=limit_[CDYN_LIMIT_ACTIVE];
        cDynJointLimit* l=list->next();
    while (l != NULL) {
        int res=l->check(t);
        if (res <= 0) { // violation
            if (res < 0) {
                list->pop(l);
                cDynTime tn=l->joint()->lastEqual(l->value(),l->error(), t, l->type());
                if (tn <= tmin) { // something not good
#if 0
                    cDynTime inc=tmax-tmin;
                    while (inc > 1e-8) { // attempt to recover by integrating at half the previous interval
                        cDynPrintf("ITERATION INTEGRATION REDUCTION JOINT %5.9f\n",inc);
                        inc *= 0.5;
                        tn=tmin+inc;
                        backup(tmin);
                        integrate(tmin,inc);
                        if (l->check(tn) >= 0) { // was able to recover
                            return(checkJointLimits(tn));
                        }
                    }
#endif
                    // can not recover. Mark joint limit as invalid
                    cDynPrintf("ITERATION FAILURE JOINT\n");
                    list->remove(l);
                    limit_[CDYN_LIMIT_INVALID]->insert(l);
#if 0
                    // restore solution to previous increment
                    inc=t-tmin;
                    backup(tmin);
                    integrate(tmin,inc);
#endif
                    // check the rest of the limits
                    return(checkJointLimits(t));
                }
                if (tn < t) { // more recent the previous limit so flush old list
                    t=tn;
                }
            }
            // add constraint to list
            contact()->jointLimit(t,l->joint(),(l->type() < 0)?CDYN_LOWER:CDYN_UPPER);
            if (flag) {
                cDynPrintf("LIMIT_J time=%-13.9f %s\n",t,l->joint()->data());
            }
            // sanity check during debugging, should be remove for final version XXXX
            assert(!(t== CDYN_TIME_UNINITIALIZED || t < tmin || t > tmax));
        }
        l=list->next(l);
    }
    return(t);
}
#endif

//---------------------------------------------------------------------------

cDynTime cDynBaseNode::checkJointLimits(const cDynTime& time)
{
    if (status_ != CDYN_ACTIVE) return(time);
    static int flag=0;
    flag=(cDynDebugLevel > 3)?1:0;
    cDynTime t=time;
    cDynTime tmin=world()->time();

    cDynJointLimitList* list=limit_[CDYN_LIMIT_ACTIVE];
        cDynJointLimit* l=list->next();
    while (l != NULL) 
    {
        int res=l->check(t);
        if (res != 1) 
        { 
            // collision has occured
            if (res != 0) 
            {	
                // a violation exists, find prior time when contact existed
                list->pop(l);
                cDynTime tu=t;
                cDynTime tl=tmin;
                cDynTime ti,tn;
                while (res != 0) 
                {
                    ti=tu-tl;
                    if (ti < 1e-8) 
                    {
                        // can not recover. Mark joint limit as invalid
                        cDynPrintf("ITERATION FAILURE JOINT\n");
                        list->remove(l);
                        limit_[CDYN_LIMIT_INVALID]->insert(l);

                        // check the rest of the limits
                        return(checkJointLimits(t));
                    }
                    tn=tl + 0.5*ti;
                    res=l->check(tn);
                    if (res == 1) tl=tn;
                    else tu=tn;
                }
                t=tu;
            }
            // add constraint to list
            contact()->jointLimit(t,l->joint(),(l->type() < 0)?CDYN_LOWER:CDYN_UPPER);
            if (flag) 
            {
                cDynPrintf("LIMIT_J time=%-13.9f %s\n",t,l->joint()->data());
            }
        }
        l=list->next(l);
    }
    return(t);
}

//---------------------------------------------------------------------------

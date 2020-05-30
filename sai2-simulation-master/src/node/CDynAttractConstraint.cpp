//===========================================================================
/*
    This file is part of the dynamics library.
    Copyright (C) 2014, Artificial Intelligence Laboratory,
    Stanford University. All rights reserved.

    \version   $MAJOR.$MINOR.$RELEASE $Rev: 1242 $
*/
//===========================================================================

//---------------------------------------------------------------------------
#include "node/CDynWorld.h"
#include "node/CDynBaseNode.h"
#include "object/CDynObject.h"
#include "node/CDynAttractConstraint.h"
//---------------------------------------------------------------------------

cDynAttractConstraint::cDynAttractConstraint(cDynObject* a, cDynVector3& pa, cDynObject* b, cDynVector3& pb, double vmax)
{
    cDynTime time=a->baseNode()->world()->time() ;  
    a_=a;
    b_=b;
    velocity_=vmax;                               
    va_.zero();
    maxImpulse_=0.0;
    maxForce_=0.0;
    
    cDynFrame fa=a->globalFrame(time);
    la_=pa;
    lb_=pb;
}

//---------------------------------------------------------------------------

void cDynAttractConstraint::check()

{
    //static int iter=0;
    //iter++;
    // use node call to indicate which objects are being considered;
    node(a_); node(b_);
    // get time interval being looked at.
    cDynTime tmin=time(0);
    cDynTime t=time(1);

    if (tmin < -1) return;

    cDynFrame fa;
    cDynFrame fb;

    fa=a_->globalFrame(t);

    cDynVector3 ga;
    cDynVector3 gb;
    ga.multiply(fa, la_);

    if (b_ != NULL) 
    {
        fb=b_->globalFrame(t);
        gb.multiply(fb, lb_);
    } 
    else 
    {
        gb=lb_;
    }

    pos_=ga;
    diff_.subtract(gb, ga);
    cDynTime dt=t-tmin;
    va_ = diff_ ;
    va_ *= 1.0/dt;
    double vmag=va_.magnitude();
    if (vmag > velocity_)
        va_ *= velocity_/vmag; 

    //cDynPrintf("VEL %5.9f (%5.9f,%5.9f,%5.9f)\n",t,va_[0],va_[1],va_[2]);

    // mark that a constraint exists between these nodes at this time
    constraint(t, a_, b_);
}

//---------------------------------------------------------------------------

void cDynAttractConstraint::define(const cDynTime& time, cDynObject* a, cDynObject* b)
{
    // add three bi-lateral constraints
    for (int i=0;i<3;i++) 
    {
        info_[i].velocity(va_[i]);
        info_[i].maxImpulse(maxImpulse_);
        info_[i].maxForce(maxForce_);
    }
    
    double maxerr = 1e-3; // doesn't matter since error is alway zero
    cDynVector3 axis;
    axis.set(1,0,0);
    //diff_.zero();
    add(axis, pos_, (fabs(diff_[0]) > maxerr)?0.0:diff_[0], maxerr, 1, &info_[0]);
    axis.set(0,1,0);
    add(axis, pos_, (fabs(diff_[1]) > maxerr)?0.0:diff_[1], maxerr, 1, &info_[1]);
    axis.set(0,0,1);
    add(axis, pos_, (fabs(diff_[2]) > maxerr)?0.0:diff_[2], maxerr, 1, &info_[2]);
}

//---------------------------------------------------------------------------
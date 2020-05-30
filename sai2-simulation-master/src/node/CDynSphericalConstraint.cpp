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
#include "node/CDynSphericalConstraint.h"
//---------------------------------------------------------------------------

cDynSphericalConstraint::cDynSphericalConstraint(cDynObject* a, cDynFrame& Fa, cDynObject* b, cDynFrame& Fb, double vmax, double wmax)
{
    a_=a;
    b_=b;
    vmax_=vmax;
    wmax_=wmax;
    va_.zero();
    maxImpulse_=0.0;
    maxForce_=0.0;

    la_=Fa;
    lb_=Fb;
}

//---------------------------------------------------------------------------

void cDynSphericalConstraint::check()
{
    //
    // static int iter=0;
    // iter++;
    // use node call to indicate which objects are being considered;
    //

    node(a_); node(b_);
    
    // get time interval being looked at.
    cDynTime tmin=time(0);
    cDynTime t=time(1);

    if (tmin < -1) return;

    cDynFrame fa;
    cDynFrame fb;

    fa=a_->globalFrame(t);
    cDynFrame ga;
    cDynFrame gb;
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
    difft_.subtract((gb.translation()), (ga.translation()));
    diffr_.angularError((gb.rotation()), (ga.rotation()));

    cDynTime dt=t-tmin;
    va_ = difft_ ;
    va_ *= 1.0/dt;
    double vmag=va_.magnitude();

    if (vmag > vmax_)
        va_ *= vmax_/vmag; 

    wa_ = diffr_;
    wa_ *= 1.0/dt;
    double wmag=wa_.magnitude();
    
    if (wmag > wmax_)
        wa_ *= wmax_/wmag;
    
    // mark that a constraint exists between these nodes at this time
    constraint(t, a_, b_);
}

//---------------------------------------------------------------------------

void cDynSphericalConstraint::define(const cDynTime& time, cDynObject* a, cDynObject* b)
{
    int i;
    
    // add three bi-lateral constraints
    for (i=0;i<3;i++) 
    {
        info_[i].velocity(va_[i]);
        info_[i].maxImpulse(maxImpulse_);
        info_[i].maxForce(maxForce_);
    }
    for (i=0;i<3;i++) 
    {
        int j=(i+1)%3;
        int k=(i+2)%3;
        info_[3+i].velocity(va_[j] + wa_[k]);
        info_[3+i].maxImpulse(1000*maxImpulse_);
        info_[3+i].maxForce(1000*maxForce_);
    }

    //	info_[3].velocity(va_[1] + wa_[2]);
    //	info_[3].velocity(va_[2] - wa_[1]);
    double maxerr = 1e-3; // doesn't matter since error is alway zero
    cDynVector3 axis;
    axis.set(1,0,0);
    //diff_.zero();
    cDynVector3 center; center=(pos_.translation());
    add(axis, center, (fabs(difft_[0]) > maxerr)?0.0:difft_[0], maxerr, 1, &info_[0]);
    axis.set(0,1,0);
    add(axis, center, (fabs(difft_[1]) > maxerr)?0.0:difft_[1], maxerr, 1, &info_[1]);
    axis.set(0,0,1);
    add(axis, center, (fabs(difft_[2]) > maxerr)?0.0:difft_[2], maxerr, 1, &info_[2]);
    cDynVector3 offset;
    offset=center; offset[0] += 1.0;

    #if 0
    axis.set(0,1,0);
    add(axis, offset, 0.0, maxerr, 1, &info_[3]);
    axis.set(0,0,1);
    add(axis, offset, 0.0, maxerr, 1, &info_[4]);
    #endif

    #if 1
    axis.set(0,1,0);
    add(axis, offset, 0.0, maxerr, 1, &info_[3]);
    offset=center; offset[1] += 1.0;
    axis.set(0,0,1);
    add(axis, offset, 0.0, maxerr, 1, &info_[4]);
    offset=center; offset[2] += 1.0;
    axis.set(1,0,0);
    add(axis, offset, 0.0, maxerr, 1, &info_[5]);
    #endif
}

//---------------------------------------------------------------------------


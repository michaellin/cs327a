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
#include "node/CDynPinConstraint.h"
//---------------------------------------------------------------------------

cDynPinConstraint::cDynPinConstraint(cDynVector3& pos, double maxerr, cDynObject* a, cDynObject* b)
{
    cDynTime time=a->baseNode()->world()->time();
    a_=a;
    b_=b;
    ga_=pos;
    gb_=pos;
    
    cDynFrame fa=a->globalFrame(time);
    la_.inversedMultiply(fa, ga_);
    if (b != NULL) 
    {
        cDynFrame fb=b->globalFrame(time);
        lb_.inversedMultiply(fb, gb_);
    }
    diff_.zero();
    maxerr_=maxerr;
    invalid_=false;
}

//---------------------------------------------------------------------------

void cDynPinConstraint::check()
{
    //static int iter=0;
    //iter++;
    // use node call to indicate which objects are being considered;
    node(a_); node(b_);
    // get time interval being looked at.
    cDynTime tmin=time(0);
    cDynTime t=time(1);

    cDynFrame fa;
    cDynFrame fb;

    double errsqr = maxerr_*maxerr_;
    double e=errsqr;

    // find time when points where less then max error away
    while (t - tmin >= 1e-8) 
    {
        fa=a_->globalFrame(t);
        ga_.multiply(fa, la_);
        if (b_ != NULL) 
        {
            fb=b_->globalFrame(t);
            gb_.multiply(fb, lb_);
        }
        diff_.subtract(gb_, ga_);
        e=diff_.dot(diff_);
        if (e < errsqr) break;
        else t=tmin+0.5*(t-tmin);
    }

    // if could not find time then mark constraint as invalid and return
    invalid_ = (t - tmin < 1e-8);
    if (invalid_) return;
    
    // else mark that a constraint exists between these nodes at this time
    constraint(t, a_, b_);
}

//---------------------------------------------------------------------------

void cDynPinConstraint::define(const cDynTime& time, cDynObject* a, cDynObject* b)
{
    // add three bi-lateral constraints
    cDynVector3 axis;
    axis.set(1,0,0);
    add(axis, ga_, diff_[0], maxerr_, 1, NULL);
    axis.set(0,1,0);
    add(axis, ga_, diff_[1], maxerr_, 1, NULL);
    axis.set(0,0,1);
    add(axis, ga_, diff_[2], maxerr_, 1, NULL);
}

//---------------------------------------------------------------------------
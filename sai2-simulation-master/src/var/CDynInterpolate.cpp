//===========================================================================
/*
    This file is part of the dynamics library.
    Copyright (C) 2014, Artificial Intelligence Laboratory,
    Stanford University. All rights reserved.

    \version   $MAJOR.$MINOR.$RELEASE $Rev: 1242 $
*/
//===========================================================================

//---------------------------------------------------------------------------
#include "matrix/CDynVector3.h"
#include "matrix/CDynVector6.h"
#include "var/CDynInterpolate.h"
//---------------------------------------------------------------------------

void cDynInterpolate::update(const double dt)
{
    t_[0].set(1,dt,dt*dt);
    t_[1].multiply(t_[0], dt*dt*dt);
    ti_=t_;

    /////////////////////////////////////////////////////////////////////////
    // Matrix:
    //      [ t^3    t^4    t^5 ]
    //  S = [ 3t^2  4t^3   5t^4 ]
    //      [ 6t   12t^2  20t^3 ]
    /////////////////////////////////////////////////////////////////////////

    cDynMatrix3 S;
    S.set(
        t_[1][0],     t_[1][1],     t_[1][2],
        3*t_[0][2],  4*t_[1][0],  5*t_[1][1],
        6*t_[0][1], 12*t_[0][2], 20*t_[1][0]
    );
    
    LU_.luDecomp(S);

    dt_=dt;
    valid_++;
}

//---------------------------------------------------------------------------

void cDynInterpolate::generate(const double dt,cDynVector6& a, const cDynVector3& ql, const cDynVector3& qu)
{
    if (dt != dt_ || !valid_) update(dt);

    a[0].set(ql[0],ql[1],0.5f*ql[2]);
    cDynVector3 y;
    y.set(
        qu[0] - a[0].dot(t_[0]),
        qu[1] - 2*a[0][2]*dt_ - a[0][1],
        -2*a[0][2]
    );
    a[1].backSub(LU_,y);
}

//---------------------------------------------------------------------------

void cDynInterpolate::solve(const double pt, const  cDynVector6& a, cDynVector3& r)
{
    if (pt != ti_[0][1]) 
    {
        ti_[0].set(1,pt,pt*pt);
        ti_[1].multiply(ti_[0], pt*pt*pt);
    }
    r.set(
        a.dot(ti_),

        5*a.elementAt(5)*ti_.elementAt(4) +
        4*a.elementAt(4)*ti_.elementAt(3) +
        3*a.elementAt(3)*ti_.elementAt(2) +
        2*a.elementAt(2)*ti_.elementAt(1) +
        a.elementAt(1),

        20*a.elementAt(5)*ti_.elementAt(3) +
        12*a.elementAt(4)*ti_.elementAt(2) +
         6*a.elementAt(3)*ti_.elementAt(1) +
         2*a.elementAt(2)
    );
}

//---------------------------------------------------------------------------
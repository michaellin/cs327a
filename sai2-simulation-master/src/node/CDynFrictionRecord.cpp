//===========================================================================
/*
    This file is part of the dynamics library.
    Copyright (C) 2014, Artificial Intelligence Laboratory,
    Stanford University. All rights reserved.

    \version   $MAJOR.$MINOR.$RELEASE $Rev: 1242 $
*/
//===========================================================================

//---------------------------------------------------------------------------
#include "node/CDynFrictionRecord.h"
//---------------------------------------------------------------------------
cDynFrictionRecord* cDynFrictionRecord::free_=NULL;
//---------------------------------------------------------------------------

cDynFrictionRecord::cDynFrictionRecord(const cDynVector3& n, cDynPrimPair* prim)
{
    rc_=1;
    
    ug_=prim->friction(CDYN_FRICTION_GRIP);
    uv_=prim->friction(CDYN_FRICTION_VISCOUS);
    ud_=prim->friction(CDYN_FRICTION_DYNAMIC);
    us_=prim->friction(CDYN_FRICTION_STATIC);
    vl_=prim->friction(CDYN_VELOCITY_LIMIT);

    if (ud_ == 0.0f && us_ == 0.0f && ug_ == 0.0f) 
    {
        friction_=false;
    } 
    else 
    {
        friction_=true;
        yt_.zero();

        //find the smallest element of n
        int i=0;
        if (cDynFabs(n[1]) < cDynFabs(n[0])) i=1;
        if (cDynFabs(n[2]) < cDynFabs(n[i])) i=2;
        yt_[i]=1;
        xt_.crossMultiply(yt_, n);
        xt_.normalize();
        yt_.crossMultiply(n, xt_);
        yt_.normalize();
    }
    
    next_ = NULL;
};

//---------------------------------------------------------------------------

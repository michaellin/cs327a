//===========================================================================
/*
    This file is part of the dynamics library.
    Copyright (C) 2014, Artificial Intelligence Laboratory,
    Stanford University. All rights reserved.

    \version   $MAJOR.$MINOR.$RELEASE $Rev: 1242 $
*/
//===========================================================================

//---------------------------------------------------------------------------
#include "node/CDynPrimPair.h"
#include "object/CDynObject.h"
#include "node/CDynConstraint.h"
//---------------------------------------------------------------------------

void cDynConstraint::constraint(cDynTime& time, cDynObject* a, cDynObject* b)
{
    if (a == NULL) return;
    cDynPrimPair* p=new cDynPrimPair(a, NULL, b, NULL);
    a->baseNode()->contact()->contact(time,p,this);
}

//---------------------------------------------------------------------------

void cDynConstraint::add(cDynVector3& n, cDynVector3& pos, double err, double maxerr, int type, cDynConstraintInfo* info)
{
    cDynContactPointType ctype=(type == 0)?CDYN_CONSTRAINT_NORMAL:CDYN_CONSTRAINT_BILATERAL;
    //cDynContactPoint* cp=
    new cDynContactPoint(contact_,n,pos,err,maxerr, ctype, NULL, info);
}

//---------------------------------------------------------------------------

void cDynConstraint::check()
{
}

//---------------------------------------------------------------------------

void cDynConstraint::define(const cDynTime& time, cDynObject* a, cDynObject* b)
{
}

//---------------------------------------------------------------------------
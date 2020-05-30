//===========================================================================
/*
    This file is part of the dynamics library.
    Copyright (C) 2014, Artificial Intelligence Laboratory,
    Stanford University. All rights reserved.

    \version   $MAJOR.$MINOR.$RELEASE $Rev: 1242 $
*/
//===========================================================================

//---------------------------------------------------------------------------
#include "object/CDynForce.h"
#include "object/CDynForceProperty.h"
//---------------------------------------------------------------------------

cDynForce::cDynForce(cDynObject* obj)
{
    properties_=NULL; obj_=obj;
}

//---------------------------------------------------------------------------

cDynForce::~cDynForce() 
{
    while (properties_)
        properties_->remove((cDynProperty **)&properties_);
}

//---------------------------------------------------------------------------

void cDynForce::add(cDynForceProperty* property) 
{
    property->insert((cDynProperty **)&properties_);
}

//---------------------------------------------------------------------------

void cDynForce::remove(cDynForceProperty* property) 
{
    property->remove((cDynProperty **)&properties_);
}

//---------------------------------------------------------------------------

void cDynForce::get(cDynVector6& force) 
{
    cDynForceRecord fr;
    fr.obj=obj_;
    fr.accum= &force;
    if (properties_) properties_->call(&fr);
}

//---------------------------------------------------------------------------
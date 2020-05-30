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
#include "object/CDynObject.h"
#include "utility/CDynError.h"
#include <assert.h>
//---------------------------------------------------------------------------

void cDynJointList::reset()
{
    while (head_ != NULL) {
        head_->remove(&head_,&tail_);
        delete head_;
    }
}

//---------------------------------------------------------------------------

cDynJointList::~cDynJointList()
{
    //remove all joint variables
    reset();
}

//---------------------------------------------------------------------------

cDynJointList::cDynJointList()
{
    obj_=NULL;
    head_=NULL;
    tail_=NULL;
}

//---------------------------------------------------------------------------

cDynJoint* cDynJointList::revolute(cDynAxis axis, char* data)
{
    assert(obj_ != NULL);
    cDynJoint* joint=new cDynJoint(obj_,CDYN_REVOLUTE,axis,data);
    joint->insert(&head_,&tail_);
    return(joint);
}

//---------------------------------------------------------------------------

cDynJoint* cDynJointList::prismatic(cDynAxis axis, char* data)
{ 
    assert(obj_ != NULL);
    cDynJoint* joint=new cDynJoint(obj_,CDYN_PRISMATIC,axis,data);
    joint->insert(&head_,&tail_);
    return(joint);
}

//---------------------------------------------------------------------------

cDynJoint* cDynJointList::spherical(char* data)
{
    assert(obj_ != NULL);
    cDynJoint *joint=new cDynJoint(obj_,CDYN_SPHERICAL,data);
    joint->insert(&head_,&tail_);
    return(joint);
}

//---------------------------------------------------------------------------

void cDynJointList::remove(cDynJoint* joint)
{
    joint->remove(&head_,&tail_);
}

//---------------------------------------------------------------------------
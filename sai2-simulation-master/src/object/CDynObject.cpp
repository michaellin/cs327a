//===========================================================================
/*
    This file is part of the dynamics library.
    Copyright (C) 2014, Artificial Intelligence Laboratory,
    Stanford University. All rights reserved.

    \version   $MAJOR.$MINOR.$RELEASE $Rev: 1242 $
*/
//===========================================================================

//---------------------------------------------------------------------------
#include "object/CDynObject.h"
//---------------------------------------------------------------------------
#include <assert.h>
#include <string.h>
//---------------------------------------------------------------------------
#include "global/CDynGlobalDefn.h"
#include "matrix/CDynMathDefn.h"
#include "node/CDynWorld.h"
#include "distance/CDynBSphere.h"
#include "node/CDynCollision.h"
#include "distance/CDynPrim.h"
#include "matrix/CDynFrame.h"
#include "node/CDynBaseNode.h"
#include "node/CDynJointLimit.h"
#include "object/CDynForceProperty.h"
//---------------------------------------------------------------------------
unsigned long cDynObject::UID_ = 0;
//---------------------------------------------------------------------------

void cDynObject::initialize(void *d)
{
    frame.reset();
    geometry.frame(&frame);
    dynamics.frame(&frame);
    joint.object(this);
    force.object(this);

    fixed_ = false;

    homeFrame_.identity();
    globalFrame_.identity();

    callnum_ = -1;

    child_=NULL;
    sibling_=NULL;
    parent_=NULL;

    collision_=NULL;

    uid_ = cDynObject::UID_++;

    baseNode_=NULL;

    data_ = d;
    controlParam_ = NULL;

    _abNode = NULL;
}

//---------------------------------------------------------------------------

cDynObject::cDynObject(void *d)
{
    initialize(d);
}

//---------------------------------------------------------------------------

cDynObject::cDynObject(char *d)
{
    initialize((void *)d);
}

//---------------------------------------------------------------------------

void cDynObject::link(cDynObject* obj)
{

    assert(obj->parent_ == NULL && obj->sibling_ == NULL); 

    cDynBaseNode* base=baseNode();
    if (base) base->beginChange();

        //adopt obj as child
        obj->parent_=this;
        obj->sibling_=child_;
        child_=obj;

        //set homeFrame of child to current frame
        obj->homeFrame_=frame.top();
    
     if (base) base->endChange();
}

//---------------------------------------------------------------------------

void cDynObject::unlink(cDynObject* obj)
{
    cDynObject* prev=NULL;
    cDynObject* p=NULL;
    for (p=child_; p != obj && p != NULL; p=p->sibling_) 
    {
        prev=p;
    } 
    assert(p == obj);

    cDynBaseNode* base=baseNode();
    if (base) base->beginChange();
        if (prev == NULL) 
        {
            child_=obj->sibling_;	
        } 
        else 
        {
            prev->sibling_=obj->sibling_;
        }
        obj->parent_=NULL;
        obj->sibling_=NULL;
        obj->baseSet(NULL,true);
    if (base) base->endChange();
}

//---------------------------------------------------------------------------

void cDynObject::unlink()
{
    assert(parent_ != NULL);
    parent_->unlink(this);
}

//---------------------------------------------------------------------------

cDynObject::~cDynObject(void)
{
    delete _abNode;

    // unlink all children
    while (child())
        unlink(child());
}

//---------------------------------------------------------------------------

void cDynObject::baseSet(cDynBaseNode* base, bool tree)
{
    cDynJoint *j;

    assert(base == NULL || baseNode() == NULL || base == baseNode());
    if (base != baseNode())
        baseNode(base);

    // add/remove joints from state space
    cDynState* state=base?base->state():NULL;
    for (j=joint.head(); j != NULL; j=j->next()) 
    {
        if (base != j->object()->baseNode()) j->object()->baseNode(base);
        cDynState* jstate=j->state()?j->state()->state():NULL;
        if (state != jstate) 
        {
            j->state(state);
        }
    }
    
    // if world is defined update state variable to last safe time
    if (state && base->world())
        state->update(base->world()->time());

    // add/remove any geometry to the collision test space
    cDynCollision* collision=base?base->collision():NULL;
    if (collision != collision_) 
    {
        if (collision_) collision_->remove(this);
        collision_=collision;
        if (collision_) collision_->add(this);
    }

    // add/remove any joint limits from base joint limit list
    cDynJointLimitList* list=base?base->limit(CDYN_LIMIT_INVALID):NULL;
    for (j=joint.head(); j != NULL; j=j->next())
        list->insert(j);

    // recurse for full tree if required
    if (tree)	
        for (cDynObject* c=child(); c != NULL; c=c->sibling())
            c->baseSet(base,tree);
}

//---------------------------------------------------------------------------

void cDynObject::collisionUpdateTree(const cDynTime& min,const cDynTime& max)
{
    if (collision_ && !isFixed())
        collision_->update(this,min,max);

    // recurse the whole tree
    for (cDynObject* c=child(); c != NULL; c=c->sibling())
        c->collisionUpdateTree(min,max);
}

//---------------------------------------------------------------------------

const cDynFrame &cDynObject::globalFrame(const cDynTime &time)
{
    if (time != baseNode()->queryTime()) 
    {
        // new time being queried, update flags to mark base as dirty

        // mark all other nodes globalFrame as dirty
        baseNode()->callinc();		

        // set the new query time
        baseNode()->queryTime(time);	
    }
    _UpdateTransformationPath(&time);
    return globalFrame_;
}

//---------------------------------------------------------------------------

void cDynObject::updateGlobalTransformation()
{
    assert(!isRoot());

    cDynFrame f = homeFrame_;
    cDynFrame fl;
    for (cDynJoint *j = joint.head(); j != NULL; j = j->next())
    {
        fl = f;
        f.multiply(fl, j->localFrame());
    }
    globalFrame_.multiply(parent_->globalFrame_, f);
}

//---------------------------------------------------------------------------

void cDynObject::_UpdateTransformationPath(const cDynTime *time)
{
    if (callnum() != baseNode()->callnum()) 
    {
        if (!isRoot() && parent_->callnum() != baseNode()->callnum())
            parent_->_UpdateTransformationPath(time);

        if (!isRoot())
        {
            for (cDynJoint *j = joint.head(); j != NULL; j = j->next())
                j->updateLocalTransformation(time);	
            updateGlobalTransformation();
        }

        callnum() = baseNode()->callnum();
    }
}

//---------------------------------------------------------------------------

#if 0
void cDynObject::extForceUpdate() // only used by DeDynamics::extForce
{
    if (forceProperty_) 
        forceProperty_->call(this);
}

//---------------------------------------------------------------------------

void cDynObject::extForce(const cDynVector6 &f)  // only used by cDynForceProperty::force()
{
    *abNode()->Fext() += f;
}
#endif

//---------------------------------------------------------------------------
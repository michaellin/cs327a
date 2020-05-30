 //===========================================================================
/*
    This file is part of the dynamics library.
    Copyright (C) 2014, Artificial Intelligence Laboratory,
    Stanford University. All rights reserved.

    \version   $MAJOR.$MINOR.$RELEASE $Rev: 1242 $
*/
//===========================================================================

//---------------------------------------------------------------------------
#include "object/CDynJoint.h"
//---------------------------------------------------------------------------
#include <assert.h>
#include <string.h>
//---------------------------------------------------------------------------
#include "utility/CDynLogger.h"
#include "object/CDynObject.h"
#include "node/CDynJointLimit.h"
#include "node/CDynBaseNode.h"
#include "var/CDynVar.h"
//---------------------------------------------------------------------------

void cDynJoint::insert(cDynJoint** head, cDynJoint** tail)
{
    prev_= *tail;
    next_=  NULL;
    if (*head == NULL) *head = this;
    if (*tail != NULL) (*tail)->next_=this;
    *tail=this;
}

//---------------------------------------------------------------------------

void cDynJoint::remove(cDynJoint** head, cDynJoint** tail)
{
    if (*head == this) *head=next_;
    if (*tail == this) *tail=prev_;
    if (next_ != NULL) next_->prev_=prev_;
    if (prev_ != NULL) prev_->next_=next_;
    obj_=NULL;
}

//---------------------------------------------------------------------------

cDynJoint::~cDynJoint()
{
    assert(obj_ != NULL);
    obj_->joint.remove(this);

    if (q_ != NULL) delete q_;
    if (sq_ != NULL) delete sq_;
    if (jcol_ != NULL) delete jcol_;
    if (hcol_ != NULL) delete hcol_;
    if (jsphere_ != NULL) delete [] jsphere_;
    if (tmp3_ != NULL) delete tmp3_;
}

//---------------------------------------------------------------------------

cDynJoint::cDynJoint(cDynObject* obj,cDynJointType type, cDynAxis axis, char* d)
{
    assert(type != CDYN_SPHERICAL);
    obj_=obj;
    type_=type;
    axis_=axis;
    data(d);

    localFrame_.identity();

    prev_=NULL;
    next_=NULL;

    sq_=NULL;
    q_ =new cDynJointVar();

    bound_[0]=bound_[1]=NULL;
    damping_=0.0f;
    inertia_=0.0f;

    jcol_=new cDynVector6;
    hcol_=new cDynVector6;
    jsphere_=NULL;

    tmp_=0.0f;
    tmp3_=NULL;

    data_ = NULL;
    controlParam_ = NULL;
}

//---------------------------------------------------------------------------

cDynJoint::cDynJoint(cDynObject* obj,cDynJointType type, char* d)
{
    assert(type == CDYN_SPHERICAL);
    obj_=obj;
    type_=type;
    //axis_=0;
    data(d);

    localFrame_.identity();

    sq_=new cDynJointSphereVar();
    q_ =NULL;

    bound_[0]=bound_[1]=NULL;
    damping_=0.0f;
    inertia_=0.0f;

    jcol_=NULL;
    hcol_=NULL;
    jsphere_=new cDynMatrix3[2];

    tmp_=0.0f;
    tmp3_=new cDynVector3;

    data_ = NULL;
    controlParam_ = NULL;
}

//---------------------------------------------------------------------------

cDynStateEntry* cDynJoint::state()
{
    return((type_ != CDYN_SPHERICAL)?q_->state():sq_->state());
}

//---------------------------------------------------------------------------

void cDynJoint::state(cDynState* s)
{
    if (type_ != CDYN_SPHERICAL)
        q_->state(s);
    else
        sq_->state(s);
}

//---------------------------------------------------------------------------

void cDynJoint::inertia(double Im) 
{
    inertia_=Im;
}

//---------------------------------------------------------------------------

void cDynJoint::damping(double b) 
{
    damping_=b;
}

//---------------------------------------------------------------------------

void cDynJoint::position(double q) 
{
    assert(q_ != NULL);

    if (q_->state())
        q_->q(q);
    else
        q_->qCurrent(q);
}

//---------------------------------------------------------------------------

void cDynJoint::velocity(double v) 
{
    assert(q_ != NULL);

    if (q_->state())
        q_->v(v);
    else
        q_->vCurrent(v);
}

//---------------------------------------------------------------------------

void cDynJoint::positionSpherical(cDynQuaternion& q) 
{
    assert(sq_ != NULL);

    if (sq_->state())
        sq_->q(q);
    else
        sq_->qCurrent(q);
}

//---------------------------------------------------------------------------

void cDynJoint::velocitySpherical(cDynVector3& w) 
{
    assert(sq_ != NULL);

    if (sq_->state())
        sq_->v(w);
    else
        sq_->vCurrent(w);
}

//---------------------------------------------------------------------------

void cDynJoint::bound(cDynBoundType type, double value)
{
    if (value == CDYN_NOBOUND) 
    {
        unbound(type);
    } 
    else 
    {
        cDynJointLimit* l=bound_[(int)type];
        if (l != NULL) 
        {
            l->value(value);
        } 
        else 
        {
            l=bound_[(int)type] = new cDynJointLimit((type == CDYN_LOWER)?-1:1,value,CDYN_MIN_JOINT_ERROR,0.0f,this);
            if (object()->baseNode()) 
            {
                    object()->baseNode()->limit(CDYN_LIMIT_INVALID)->insert(l);
            }
        }
    }
}

//---------------------------------------------------------------------------

void cDynJoint::unbound(const cDynBoundType type)
{
    cDynJointLimit* l=bound_[(int)type];
    if (l != NULL) 
    {
        int state=l->state();
        if (state >= 0)
            object()->baseNode()->limit(state)->remove(l);

        delete l;
        bound_[(int)type]=NULL;
    }
}

//---------------------------------------------------------------------------

void cDynJoint::invalid(const cDynBoundType type)
{
    cDynJointLimit* l=bound_[(int)type];
    assert(l != NULL);
    if (l->state() != CDYN_LIMIT_INVALID) 
    {
        if (l->state() >= 0)
            object()->baseNode()->limit(l->state())->remove(l);

        object()->baseNode()->limit(CDYN_LIMIT_INVALID)->insert(l);
        cDynPrintf("Joint %s %c marked as valid\n",data(),(type == CDYN_LOWER)?'<':'>');
    }
}

//---------------------------------------------------------------------------

void cDynJoint::error(const cDynBoundType type,const double value)
{
    if (bound_[(int)type] != NULL) bound_[(int)type]->error(value);
}

//---------------------------------------------------------------------------

void cDynJoint::epsilon(const cDynBoundType type,const double value)
{
    if (bound_[(int)type] != NULL) bound_[(int)type]->epsilon(value);
}

//---------------------------------------------------------------------------

double cDynJoint::bound(cDynBoundType type) const
{
    if (bound_[(int)type] == NULL) return(CDYN_NOBOUND);
    else return(bound_[(int)type]->value());
}

//---------------------------------------------------------------------------

double cDynJoint::error(cDynBoundType type) const
{
    if (bound_[(int)type] == NULL) return(CDYN_NOBOUND);
    else return(bound_[(int)type]->error());
}

//---------------------------------------------------------------------------

double cDynJoint::epsilon(const cDynBoundType type) const
{
    if (bound_[(int)type] == NULL) return(CDYN_NOBOUND);
    else return(bound_[(int)type]->epsilon());
}

//---------------------------------------------------------------------------

void cDynJoint::updateLocalTransformation(const cDynTime *time)
{
    switch(type())
    {	
        case CDYN_SPHERICAL:
        localFrame_.set(time ? positionSphericalAt(*time) : sq());
        break;
    
        case CDYN_PRISMATIC:
        // trans_local = transHome + RHome * q
        localFrame_.translation(axis()) = time ? positionAt(*time) : q();
        break;
    
        case CDYN_REVOLUTE:
#if 1
            cDynQuaternion tmpQ;
            tmpQ.set(axis(), time ? positionAt(*time) : q());
            localFrame_.set(tmpQ);
#else
        // R_local = RHome * R_theta
        cSetQ4S2(localFrame_.rotation(), axis(), q());
#endif
        break;
    }
}

//---------------------------------------------------------------------------

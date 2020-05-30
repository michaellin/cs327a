//===========================================================================
/*
    This file is part of the dynamics library.
    Copyright (C) 2014, Artificial Intelligence Laboratory,
    Stanford University. All rights reserved.

    \version   $MAJOR.$MINOR.$RELEASE $Rev: 1242 $
*/
//===========================================================================

//---------------------------------------------------------------------------
#include "dynamics/CDynDynamics.h"
#include "object/CDynObject.h"
#include "object/CDynJoint.h"
#include "node/CDynBaseNode.h"
#include "utility/CDynLogger.h"
//---------------------------------------------------------------------------
#include "dynamics/CDynABDynamics.h"
#include "dynamics/CDynABJoint.h"
//---------------------------------------------------------------------------
#define VALUE_ALL	0
#define VALUE_Q		1
#define VALUE_DQ	2
#define VALUE_DDQ	3
#define VALUE_TAU	4
#define VALUE_JG	5
#define VALUE_DAMPING 6
#define VALUE_DQ_TAU 7
//---------------------------------------------------------------------------

void _cDynSetValues(cDynObject* obj, int flag, int rotate)
{
    {
        int i = 0;
        for (cDynJoint* j = obj->joint.head(); j != NULL; j = j->next())
        {
            if (j->type() == CDYN_SPHERICAL) //change + optimize
            {
                cDynABJointSpherical* abJoint = (cDynABJointSpherical*)obj->abNode()->abJoint(i++);

                if (flag == VALUE_ALL || flag == VALUE_Q)
                    abJoint->_Q = j->sq();
                if (flag == VALUE_ALL || flag == VALUE_DQ || flag == VALUE_DQ_TAU)
                {
                    abJoint->_dQ = j->sv();
                    if (rotate)
                    {
                        cDynVector3 tmpV;
                        j->rotateHome2Local(tmpV, j->sv());
                        abJoint->_dQ = tmpV;
                    }
                }
                if (flag == VALUE_ALL || flag == VALUE_TAU || flag == VALUE_DQ_TAU)
                {
                    if (!rotate)
                    {
                        abJoint->_ddQ = j->torqueSpherical();
                        cDynVector3 tmpV;
                        tmpV.zero();
                        j->torqueSpherical(tmpV);
                    }
                    abJoint->_Tau = j->storqueTotal();
                }
                if (flag == VALUE_ALL || flag == VALUE_DAMPING)
                {
                    abJoint->_damping = j->damping();
                    abJoint->_motorInertia = j->inertia();
                }
            }
            else
            {
                cDynABJointDOF1* abJoint = (cDynABJointDOF1*)obj->abNode()->abJoint(i++);

                if (flag == VALUE_ALL || flag == VALUE_Q)
                    abJoint->_Q = j->q();
                if (flag == VALUE_ALL || flag == VALUE_DQ || flag == VALUE_DQ_TAU)
                    abJoint->_dQ = j->v();
                if (flag == VALUE_ALL || flag == VALUE_TAU || flag == VALUE_DQ_TAU)
                {
                    if (!rotate)
                    {
                        abJoint->_ddQ = j->torque();
                        j->torque(0);
                    }
                    abJoint->_Tau = j->torqueTotal();
                }
                if (flag == VALUE_ALL || flag == VALUE_DAMPING)
                {
                    abJoint->_damping = j->damping();
                    abJoint->_motorInertia = j->inertia();
                }
            }
        }
    }
    for (cDynObject *o = obj->child(); o != NULL; o = o->sibling())
        _cDynSetValues(o, flag, rotate);
}

//---------------------------------------------------------------------------

void _cDynGetValues(cDynObject *obj, int flag, int rotate)
{
    {
        int i = 0;
        for (cDynJoint* j = obj->joint.head(); j != NULL; j = j->next())
        {

            if (j->type() == CDYN_SPHERICAL) //change + optimize
            {
                cDynABJointSpherical* abJoint = (cDynABJointSpherical*)obj->abNode()->abJoint(i++);

                if (flag == VALUE_ALL || flag == VALUE_DDQ)
                {
                    cDynVector3 tmpV = abJoint->_ddQ;
                    if (rotate)
                    {
                        j->rotateLocal2Home(tmpV, abJoint->_ddQ);
                    }
                    j->sa(tmpV);
                }
                if (flag == VALUE_DQ)
                {
                    assert(rotate);
                    cDynVector3 tmpV = abJoint->_ddQ;
                    j->rotateLocal2Home(tmpV, abJoint->_ddQ);
                    tmpV += j->sv();
                    j->sv(tmpV);
                }
                if (flag == VALUE_TAU)
                {
                    assert(!rotate);
                    abJoint->_Tau.zero();
                    j->torqueSpherical(abJoint->_Tau);
                    abJoint->_Tau.zero();
                }
                if (flag == VALUE_ALL || flag == VALUE_JG)
                {
                    j->jsphere()[0] = abJoint->_Jg[0];
                    j->jsphere()[1] = abJoint->_Jg[1];
                }
            }
            else
            {
                cDynABJointDOF1* abJoint = (cDynABJointDOF1*)obj->abNode()->abJoint(i++);

                if (flag == VALUE_ALL || flag == VALUE_DDQ)
                    j->a(abJoint->_ddQ);

                if (flag == VALUE_DQ)
                {
                    assert(rotate);
                    j->v(j->v() + abJoint->_ddQ);
                }
                
                if (flag == VALUE_TAU)
                {
                    assert(!rotate);
                    abJoint->_Tau = 0;
                    j->torque(abJoint->_Tau);
                    abJoint->_Tau = 0;
                }
                
                if (flag == VALUE_ALL || flag == VALUE_JG)
                    j->jcol() = abJoint->_Jg;
            }
        }
    }
    for (cDynObject *o = obj->child(); o != NULL; o = o->sibling())
        _cDynGetValues(o, flag, rotate);
}

//---------------------------------------------------------------------------

void _cDynNodeUpdateStateConfigTree2(cDynObject *obj, long callnum)
{
    if (!obj->isRoot())
    {
        cDynJoint *j = obj->joint.head();
        if (j)
        {
            j->updateLocalTransformation(NULL);
            cDynFrame f;
            f.multiply(obj->homeFrame(),j->localFrame());
            j->localX().set(f);
            while (j->next())
            {
                j = j->next();
                j->updateLocalTransformation(NULL);
                j->localX().set(j->localFrame());
            }
        }
        obj->updateGlobalTransformation();
        obj->callnum() = callnum;
    }

    for (cDynObject *o = obj->child(); o != NULL; o = o->sibling())
        _cDynNodeUpdateStateConfigTree2(o, callnum);
}

//---------------------------------------------------------------------------

void _cDynNodeInitializeNodeTree2(cDynObject *obj)
{
    if (obj->isRoot()) 
    {
        obj->fixed(true); // NEED TO FIX DeDynamics shouldn't modify items in cDynNode
    } 
    else 
    {
        if (obj->parent()->isFixed() && obj->joint.head() == NULL)
            obj->fixed(true);
        else
            obj->fixed(false);
    }
    cDynJoint *j = obj->joint.head();
    if (j)
    {
        obj->localX().identity();
        while (j)
        {
            if (j->type() == CDYN_PRISMATIC)
                if (j->axis() == 0)
                    j->v(0.0);
            if (j->type() == CDYN_SPHERICAL)
            {
                cDynVector3 tmpV;
                tmpV.set(0,0,0);
                j->sv(tmpV);
            }
            j = j->next();
        }
    }
    else
    {
        obj->localX().set(obj->homeFrame());
    }
    for (cDynObject *o = obj->child(); o != NULL; o = o->sibling())
        _cDynNodeInitializeNodeTree2(o);
}

//---------------------------------------------------------------------------

void _cDynNodeInitializeJointTree(cDynObject *obj)
{
    if (obj->isRoot())
        obj->abNode(new cDynABNodeRoot);
    else
    {
        cDynABNode* nnode = NULL;
        if (obj->isRoot())
            nnode = new cDynABNodeRoot;
        else if (obj->joint.head() == obj->joint.tail())	
            nnode = new cDynABNodeNOJ1;
        else
            nnode = new cDynABNodeNOJn;

        obj->abNode(nnode);
        int i = 0;
        cDynJoint* j;
        for (j = obj->joint.head(); j != NULL; j = j->next())
            i++;
        nnode->noj(i);

        cDynABJoint *joint=NULL;

        i = 0;
        for (j = obj->joint.head(); j != NULL; j = j->next())
        {
            if (j->type() == CDYN_SPHERICAL)
                joint = new cDynABJointSpherical;
            else if (j->type() == CDYN_PRISMATIC)
                joint = new cDynABJointPrismatic(j->axis());
            else if (j->type() == CDYN_REVOLUTE)
                joint = new cDynABJointRevolute(j->axis());
                
            nnode->abJoint(joint, i);
            i++;
        }
        if (i == 0)
            obj->abNode()->abJoint(new cDynABJointFixed);
    }

    for (cDynObject *o = obj->child(); o != NULL; o = o->sibling())
        _cDynNodeInitializeJointTree(o);
}

//---------------------------------------------------------------------------

void cDynamicsInitialize2(cDynObject *root, long callnum)
{
    _cDynNodeInitializeJointTree(root);
    _cDynNodeInitializeNodeTree2(root);
    _cDynNodeUpdateStateConfigTree2(root, callnum);
    _cDynSetValues(root, VALUE_ALL, true);
    cDynABDynamics::updateLocalXTreeOut(root);
    cDynABDynamics::resetInertiaTreeOut(root);
    cDynABDynamics::resetFlagTreeOut(root);
}

//---------------------------------------------------------------------------

void cDynamicsInvDynamics(cDynObject *root, cDynVector3 *gravity, const int compensateGravity)
{
    _cDynSetValues(root, VALUE_TAU, false);
//	cDynABDynamics::updateLocalXTreeOut(root);
    cDynVector3 g;
#ifndef BORLAND
    g.inversedMultiply(root->globalFrame().rotation(), *gravity);
#else
    g.inversedMultiply(*(cDynQuaternion*)(root->globalFrame().rotation()), *gravity);
#endif


    cDynABDynamics::inverseDynamics(root, &g, compensateGravity);
    _cDynGetValues(root, VALUE_TAU, false);
}

//---------------------------------------------------------------------------

void cDynamicsFwdDynamics2(cDynObject *root, cDynVector3 *gravity, long callnum)
{
    _cDynNodeUpdateStateConfigTree2(root, callnum);
    _cDynSetValues(root, VALUE_ALL, true);
    cDynABDynamics::updateLocalXTreeOut(root);
    cDynVector3 g;
    g.inversedMultiply((root->globalFrame().rotation()), *gravity);
    cDynABDynamics::forwardDynamics(root, &g);
//	cDynABDynamics::forwardDynamicsControl(root);
    _cDynGetValues(root, VALUE_DDQ, true);
}

//---------------------------------------------------------------------------

void cDynamicsInitFwdDynamicsConfig2(cDynObject *root, long callnum)
{
    _cDynNodeUpdateStateConfigTree2(root, callnum);
    _cDynSetValues(root, VALUE_ALL, true);
    cDynABDynamics::updateLocalXTreeOut(root);
    cDynABDynamics::forwardDynamicsConfigInit(root);
}

//---------------------------------------------------------------------------

void _cDynamicsFwdDynamicsIn2(cDynObject *contactNode, double (*fwdIn)(cDynJoint *j), const cDynVector3& (*fwdInSphere)(cDynJoint *j))
{
    for (cDynJoint* j = contactNode->joint.tail(); j != NULL; j = j->prev())
    {
        if (j->type() == CDYN_SPHERICAL)
            fwdInSphere(j);
        else
            fwdIn(j);
    }
    cDynObject *o = contactNode->parent();
    if (o)
        _cDynamicsFwdDynamicsIn2(o, fwdIn, fwdInSphere);
}

//---------------------------------------------------------------------------

void cDynamicsFwdDynamicsIn2(cDynObject *contactNode, double (*fwdIn)(cDynJoint *j), const cDynVector3& (*fwdInSphere)(cDynJoint *j))
{
//	cDynamicsInitFwdDynamicsConfig2(root, root->baseNode()->callnum());

    _cDynamicsFwdDynamicsIn2(contactNode, fwdIn, fwdInSphere);
}

//---------------------------------------------------------------------------

void cDynamicsFwdDynamicsVel2(cDynObject *root)
{	
//	cDynamicsInitFwdDynamicsConfig2(root, root->baseNode()->callnum());

    _cDynSetValues(root, VALUE_DQ_TAU, true);

    // FIX: necessary?
//	_cDynNodeUpdateStateConfigTree2(root, root->baseNode()->callnum());
//	_cDynSetValues(root, VALUE_ALL, true);
//	cDynABDynamics::updateLocalXTreeOut(root);

    cDynVector3 g;
    g.inversedMultiply((root->globalFrame().rotation()), root->baseNode()->gravity);
    cDynABDynamics::forwardDynamics(root, &g);
    _cDynGetValues(root, VALUE_DDQ, false);
}

//---------------------------------------------------------------------------

void _cDynSetValuesConfig(cDynObject *obj, double (*fwdIn)(cDynJoint *j), const cDynVector3& (*fwdInSphere)(cDynJoint *j))
{
    {
        int i = 0;	
        for (cDynJoint* j = obj->joint.head(); j != NULL; j = j->next())
        {
            if (j->type() == CDYN_SPHERICAL)
            {
                cDynABJointSpherical* abJoint = (cDynABJointSpherical*)obj->abNode()->abJoint(i++);
                abJoint->_Tau = fwdInSphere(j);
            }
            else
            {
                cDynABJointDOF1* abJoint = (cDynABJointDOF1*)obj->abNode()->abJoint(i++);
                abJoint->_Tau = fwdIn(j);
            }
        }
    }
    cDynObject *o = obj->parent();
    if (!o->isRoot())
        _cDynSetValuesConfig(o, fwdIn, fwdInSphere);
}

//---------------------------------------------------------------------------

void _cDynGetValuesConfig(cDynObject *obj, void (*fwdOut)(cDynJoint* j, const double x), void (*fwdOutSphere)(cDynJoint* j, const cDynVector3& x), int rotate)
{
    {
        int i = 0;
        for (cDynJoint* j = obj->joint.head(); j != NULL; j = j->next())
        {
            if (j->type() == CDYN_SPHERICAL) //change + optimize
            {
                cDynABJointSpherical* abJoint = (cDynABJointSpherical*)obj->abNode()->abJoint(i++);

                cDynVector3 tmpV;
                if (rotate)
                    j->rotateLocal2Home(tmpV, abJoint->_ddQ);
                else
                    tmpV = abJoint->_ddQ;

                fwdOutSphere(j, tmpV);
            }
            else
            {
                cDynABJointDOF1* abJoint = (cDynABJointDOF1*)obj->abNode()->abJoint(i++);
                fwdOut(j, abJoint->_ddQ);
            }
        }
    }
    for (cDynObject *o = obj->child(); o != NULL; o = o->sibling())
        _cDynGetValuesConfig(o, fwdOut, fwdOutSphere, rotate);
}

//---------------------------------------------------------------------------

void cDynamicsFwdDynamicsConfig2(cDynObject *root,cDynObject *contact,
                                 double (*fwdIn)(cDynJoint *j), const cDynVector3& (*fwdInSphere)(cDynJoint *j),
                                 void (*fwdOut)(cDynJoint* j, const double x), void (*fwdOutSphere)(cDynJoint* j, const cDynVector3& x), int rotate)
{
    assert(!contact->isRoot());
//	cDynamicsInitFwdDynamicsConfig2(root, root->baseNode()->callnum());

    _cDynSetValuesConfig(contact, fwdIn, fwdInSphere);
    cDynABDynamics::forwardDynamicsConfig(root, contact);	
    _cDynGetValuesConfig(root, fwdOut, fwdOutSphere, rotate);
}

//---------------------------------------------------------------------------

void cDynamicsGlobalJacobian2(cDynObject *obj)
{
    cDynABDynamics::globalJacobian(obj);

    _cDynGetValues(obj, VALUE_JG, false);
}

//---------------------------------------------------------------------------

void cDynamicsTorqueExternalZeroInwardPath2(cDynObject *obj)
{
    for (cDynJoint* j = obj->joint.tail(); j != NULL; j = j->prev())
    {
        if (j->type() == CDYN_SPHERICAL)
            j->storqueExternalZero();
        else
            j->torqueExternalZero();
    }

    cDynObject *o = obj->parent();
    if (o != NULL)
        cDynamicsTorqueExternalZeroInwardPath2(o);
}

//---------------------------------------------------------------------------

double cDynamicsPotentialEnergy2(cDynObject *rootNode,const cDynVector3 *gravity)
{
    return cDynABDynamics::potentialEnergy(rootNode, gravity);
}

//---------------------------------------------------------------------------

double cDynamicsKineticEnergy2(cDynObject *rootNode)
{
    return cDynABDynamics::kineticEnergy(rootNode);
}

//---------------------------------------------------------------------------

#ifdef CDYN_DEBUG
void cDynamicsDisplay2(cDynObject *rootNode)
{
//	assert(false);
}
#endif // CDYN_DEBUG

void cDynamicsResetMass2(cDynObject *obj,const double mass)
{
    obj->dynamics.scale(mass);
    cDynABDynamics::resetInertia(obj);
}

//---------------------------------------------------------------------------

void cDynamicsImpulse2(cDynObject *rootNode,cDynObject *contact,const cDynVector3 *contactPoint,const cDynVector3 *impulseVector)
{
    assert(!contact->isRoot());

    cDynABDynamics::forwardDynamicsImpulse(rootNode, contact, contactPoint, impulseVector);
    _cDynGetValues(rootNode, VALUE_DQ, true);
}

//---------------------------------------------------------------------------

void cDynamicsInitialize(cDynObject *rootNode, long callnum)
{
    cDynamicsInitialize2(rootNode, callnum);
}

//---------------------------------------------------------------------------

void cDynamicsFwdDynamics(cDynObject *rootNode, cDynVector3 *gravity, long callnum)
{
    cDynamicsFwdDynamics2(rootNode, gravity, callnum);
}

//---------------------------------------------------------------------------

void cDynamicsInitFwdDynamicsConfig(cDynObject *rootNode, long callnum)
{
    cDynamicsInitFwdDynamicsConfig2(rootNode, callnum);
}

//---------------------------------------------------------------------------

void cDynamicsFwdDynamicsConfig(cDynObject *rootNode,cDynObject *contactNode,
                                 double (*fwdIn)(cDynJoint *j), const cDynVector3& (*fwdInSphere)(cDynJoint *j),
                                 void (*fwdOut)(cDynJoint* j, const double x), void (*fwdOutSphere)(cDynJoint* j, const cDynVector3& x), int rotate)
{
    cDynamicsFwdDynamicsConfig2(rootNode, contactNode,
                                 fwdIn, fwdInSphere,
                                 fwdOut, fwdOutSphere, rotate);
}

//---------------------------------------------------------------------------

/*
inline void cDynamicsFwdDynamicsImpulse(cDynObject *rootNode,cDynObject *impulseNode)
{
    cDynamicsFwdDynamicsImpulse2(rootNode,impulseNode);
}

*/
//---------------------------------------------------------------------------

void cDynamicsFwdDynamicsVel(cDynObject *rootNode)
{
    cDynamicsFwdDynamicsVel2(rootNode);
}

//---------------------------------------------------------------------------

void cDynamicsImpulse(cDynObject *rootNode,cDynObject *contactNode,const cDynVector3 *contactPoint,const cDynVector3 *impulseVector)
{
    cDynamicsImpulse2(rootNode,contactNode,contactPoint,impulseVector);
}

//---------------------------------------------------------------------------

void cDynamicsFwdDynamicsIn(cDynObject *contactNode, double (*fwdIn)(cDynJoint *j), const cDynVector3& (*fwdInSphere)(cDynJoint *j))
{
    cDynamicsFwdDynamicsIn2(contactNode, fwdIn, fwdInSphere);
}

//---------------------------------------------------------------------------

double cDynamicsPotentialEnergy(cDynObject *rootNode,const cDynVector3 *gravity)
{ 
    return cDynamicsPotentialEnergy2(rootNode,gravity);
}

//---------------------------------------------------------------------------

double cDynamicsKineticEnergy(cDynObject *rootNode)
{
    return cDynamicsKineticEnergy2(rootNode);
}

//---------------------------------------------------------------------------

void cDynamicsGlobalJacobian(cDynObject *rootNode)
{
    cDynamicsGlobalJacobian2(rootNode);
}

//---------------------------------------------------------------------------

void cDynamicsTorqueExternalZeroInwardPath(cDynObject *node)
{
    cDynamicsTorqueExternalZeroInwardPath2(node);
}

//---------------------------------------------------------------------------

void cDynamicsResetMass(cDynObject *node,const double mass)
{
    cDynamicsResetMass2(node,mass);
}

//---------------------------------------------------------------------------
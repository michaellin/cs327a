//===========================================================================
/*
    This file is part of the dynamics library.
    Copyright (C) 2014, Artificial Intelligence Laboratory,
    Stanford University. All rights reserved.

    \version   $MAJOR.$MINOR.$RELEASE $Rev: 1242 $
*/
//===========================================================================

//---------------------------------------------------------------------------
#include "dynamics/CDynABDynamics.h"
//---------------------------------------------------------------------------
#include "matrix/CDynVector3.h"
#include "matrix/CDynVector6.h"
#include "matrix/CDynMatrix6.h"
#include "dynamics/CDynABNode.h"
//---------------------------------------------------------------------------

class cDynABDynamicsData
{
    //-----------------------------------------------------------------------
    // PUBLIC MEMBERS:
    //-----------------------------------------------------------------------
public:

    cDynMatrix6 Ia;
    cDynVector3 WxV;
    cDynVector3 g;
};

//---------------------------------------------------------------------------

class cDynABDynamicsData2
{
    //-----------------------------------------------------------------------
    // PUBLIC MEMBERS:
    //-----------------------------------------------------------------------
public:

    cDynVector3 WxV;
    cDynVector3 g;
};

//---------------------------------------------------------------------------

void cDynABDynamics::updateLocalXTreeOut(cDynDNode* root)
{
    root->abNode()->updateLocalX(root->homeFrame());

    for (cDynDNode* n = root->child(); n != NULL; n = n->sibling())
        updateLocalXTreeOut(n);
}
//---------------------------------------------------------------------------

void cDynABDynamics::resetInertia(cDynDNode* obj)
{
    assert(obj->mass());

    obj->abNode()->inertia(obj->mass()->mass(), obj->mass()->center(), obj->mass()->inertia());
}

//---------------------------------------------------------------------------

void cDynABDynamics::resetInertiaTreeOut(cDynDNode* root)
{
    if (root->mass())
        resetInertia(root);

    for (cDynDNode* n = root->child(); n != NULL; n = n->sibling())
        resetInertiaTreeOut(n);
}

//---------------------------------------------------------------------------

void cDynABDynamics::resetFlagTreeOut(cDynDNode* root)
{
    root->abNode()->flag(0);

    for (cDynDNode* n = root->child(); n != NULL; n = n->sibling())
        resetFlagTreeOut(n);
}

//---------------------------------------------------------------------------

void cDynABDynamics::forwardDynamics(cDynDNode* root, const cDynVector3* gravity)
{
    cDynABDynamicsData data;
    data.g = *gravity;

    _forwardDynamicsOutIn(root, &data, NULL, NULL);

    _accelerationTreeOut(root, NULL);
}

//---------------------------------------------------------------------------

void cDynABDynamics::inverseDynamics(cDynDNode* root, const cDynVector3* gravity, const int compensateGravity)
{
    cDynABDynamicsData2 data;
    data.g = *gravity;

    _inverseDynamicsOutIn(root, &data, NULL, NULL, NULL, compensateGravity);
}

//---------------------------------------------------------------------------

void cDynABDynamics::forwardDynamicsConfigInit(cDynDNode* root)
{
    _forwardDynamicsConfigInitOutIn(root, NULL);
}

//---------------------------------------------------------------------------

void cDynABDynamics::forwardDynamicsConfig(cDynDNode* root, cDynDNode* contact)
{
    assert(!contact->isRoot());

    contact->abNode()->biasForceConfigInit();
    contact->abNode()->abBiasForceConfigInit();
    _articulatedBiasForceConfigPathIn(contact);

    _velocityDeltaTreeOut(root, NULL);
}

//---------------------------------------------------------------------------

void cDynABDynamics::forwardDynamicsImpulse(cDynDNode* root, cDynDNode* contact, const cDynVector3* point, const cDynVector3* impulse)
{
    assert(!contact->isRoot());

    contact->abNode()->impulseInit(*point, *impulse);
    contact->abNode()->abImpulseInit();
    _articulatedImpulsePathIn(contact);

    _velocityDeltaTreeOut(root, NULL);
}

//---------------------------------------------------------------------------

void cDynABDynamics::globalJacobian(cDynDNode* root)
{
    root->abNode()->globalJacobian(root->globalFrame());

    for (cDynDNode* n = root->child(); n != NULL; n = n->sibling())
        globalJacobian(n);
}

//---------------------------------------------------------------------------

void cDynABDynamics::_articulatedImpulsePathIn(cDynDNode* contact)
{
    cDynDNode *n = contact->parent();

    contact->abNode()->abImpulse(*(n->abNode()->Pa()), !n->isRoot());

    if (!n->isRoot())
    {
        n->abNode()->abImpulseInit();
        _articulatedImpulsePathIn(n);
    }
}

//---------------------------------------------------------------------------

void cDynABDynamics::_articulatedBiasForceConfigPathIn(cDynDNode* contact)
{
    cDynDNode *n = contact->parent();
    
    contact->abNode()->abBiasForceConfig(*(n->abNode()->Pa()), !n->isRoot());

    if (!n->isRoot())
    {
        n->abNode()->abBiasForceConfigInit();
        _articulatedBiasForceConfigPathIn(n);
    }
}

//---------------------------------------------------------------------------

void cDynABDynamics::_forwardDynamicsConfigInitOutIn(cDynDNode* root, cDynMatrix6* Iah)
{
    cDynMatrix6 Ia;

    root->abNode()->abInertiaInit(Ia);

    for (cDynDNode* n = root->child(); n != NULL; n = n->sibling())
        _forwardDynamicsConfigInitOutIn(n, &Ia);

    root->abNode()->abInertiaDependConfig(*Iah, Ia, !root->isParentRoot());
}

//---------------------------------------------------------------------------

void cDynABDynamics::_forwardDynamicsOutIn(cDynDNode* root, cDynABDynamicsData* datah, cDynVector6* Pah, const cDynVector6 *Vh)
{
    cDynABDynamicsData data;
    cDynVector6* Pa = root->abNode()->Pa();
    cDynVector6* V = root->abNode()->V();
    cDynVector6 Fext;

    root->abNode()->abInertiaInit(data.Ia);

    root->abNode()->velocity(*V, data.WxV, *Vh, datah->WxV);

    root->abNode()->biasForce(*Pa, *V, data.WxV);

    root->abNode()->gravityForce(Fext, data.g, datah->g);

    root->force.get(Fext);

    root->abNode()->externalForce(*Pa, Fext);

    for (cDynDNode* n = root->child(); n != NULL; n = n->sibling())
        _forwardDynamicsOutIn(n, &data, Pa, V);

    root->abNode()->abInertiaDepend(datah->Ia, *Pah, data.Ia, !root->isParentRoot());
}

//---------------------------------------------------------------------------

void cDynABDynamics::_inverseDynamicsOutIn(cDynDNode* root, cDynABDynamicsData2* datah, cDynVector6* Fh, const cDynVector6* Vh, const cDynVector6* Ah, const int compensateGravity)
{
    cDynABDynamicsData2 data;
    cDynVector6* F = root->abNode()->Pa();
    cDynVector6* A = root->abNode()->A();
    cDynVector6* V = root->abNode()->V();
    cDynVector6 P;
    cDynVector6 Fext;

    F->zero();
    
    if (compensateGravity == CDYN_ABDYNAMICS_CONTROL_FLOAT)
    {
        root->abNode()->gravityForce(*F, data.g, datah->g);
        F->negate(*F);
    }
    else
    {
        root->abNode()->velocity(*V, data.WxV, *Vh, datah->WxV);

        root->abNode()->biasForce(P, *V, data.WxV);
        
        root->abNode()->gravityForce(Fext, data.g, datah->g);
        
        root->abNode()->externalForce(P, Fext);

        root->abNode()->accelerationOnly(*A, *Ah);

        root->abNode()->netForce(*F, *A, P);

        if (!compensateGravity)
        {
            P.negate(Fext);
            root->abNode()->externalForce(*F, P);
        }
    }

    for (cDynDNode* n = root->child(); n != NULL; n = n->sibling())
        _inverseDynamicsOutIn(n, &data, F, V, A, compensateGravity);

    root->abNode()->force(*Fh, !root->isParentRoot());
}

//---------------------------------------------------------------------------

void cDynABDynamics::_accelerationTreeOut(cDynDNode* root, const cDynVector6* Ah)
{
    cDynVector6* A = root->abNode()->A();

    root->abNode()->acceleration(*A, *Ah);

    for (cDynDNode* n = root->child(); n != NULL; n = n->sibling())
        _accelerationTreeOut(n, A);
}

//---------------------------------------------------------------------------

void cDynABDynamics::_velocityDeltaTreeOut(cDynDNode* root, const cDynVector6* dVh)
{
    cDynVector6 dV;

    root->abNode()->velocityDelta(dV, *dVh);
    
    for (cDynDNode* n = root->child(); n != NULL; n = n->sibling())
        _velocityDeltaTreeOut(n, &dV);
}

//---------------------------------------------------------------------------

double cDynABDynamics::potentialEnergy(cDynDNode* root, const cDynVector3* gh)
{
    double E = 0;
    cDynVector3 g;

    if (root->mass())
        E = root->abNode()->potentialEnergy(g, *gh, root->globalFrame(), root->mass()->mass(), root->mass()->center());

    for (cDynDNode* n = root->child(); n != NULL; n = n->sibling())
        E += potentialEnergy(n, &g);

    return E;
}

//---------------------------------------------------------------------------

double cDynABDynamics::kineticEnergy(cDynDNode* root, const cDynVector6* Vh)
{
    double E = 0;
    cDynVector6* V = root->abNode()->V();

    if (root->mass())
        E = root->abNode()->kineticEnergy(*V, *Vh);

    for (cDynDNode* n = root->child(); n != NULL; n = n->sibling())
        E += kineticEnergy(n, V);

    return E;
}

//---------------------------------------------------------------------------


//===========================================================================
/*
    This file is part of the dynamics library.
    Copyright (C) 2014, Artificial Intelligence Laboratory,
    Stanford University. All rights reserved.

    \version   $MAJOR.$MINOR.$RELEASE $Rev: 1242 $
*/
//===========================================================================

//---------------------------------------------------------------------------
#include "global/CDynGlobalDefn.h"
#include "dynamics/CDynABNode.h"
#include "dynamics/CDynABJoint.h"
//---------------------------------------------------------------------------
#include <assert.h>
//---------------------------------------------------------------------------

void cDynABNodeNOJ::abInertiaInit(cDynMatrix6& Ia)
{ 
    Ia = *I();
}

//---------------------------------------------------------------------------

void cDynABNodeNOJ::impulseInit(const cDynVector3& point, const cDynVector3& impulse)
{
    cDynVector6& Ya = *Pa();

    // Ya = -Xc Yc
    //    = -[R , 0; dxR , R] [ y, 0]
    //    = [-Ry ; -dxRy]
    Ya[0].negate(impulse);
    Ya[1].crossMultiply(point, Ya[0]);
}

//---------------------------------------------------------------------------

void cDynABNodeNOJ::biasForceConfigInit()
{
    Pa()->zero();
}

//---------------------------------------------------------------------------

void cDynABNodeNOJ::abBiasForceConfigInit()
{
    _flag = true;
}

//---------------------------------------------------------------------------

void cDynABNodeNOJ::netForce(cDynVector6& F, const cDynVector6& A, const cDynVector6& P)
{
    F.multiply(*I(), A);
    F += P;
}

//---------------------------------------------------------------------------
//
// Pa -= Fext
//
void cDynABNodeNOJ::externalForce(cDynVector6& Pa, const cDynVector6& Fext)
{
    Pa -= Fext;
}

//---------------------------------------------------------------------------
//
// Iah = Ih + sum [ Li Iai Lti ]
//
void cDynABNodeNOJ::_abInertia(cDynMatrix6& Iah, const cDynMatrix6& L, const cDynMatrix6& Ia)
{
    cDynMatrix6 tmpM6;
    tmpM6.similarityXform(L, Ia);
    Iah += tmpM6;
}

//---------------------------------------------------------------------------
//
// Pah = Ph - Fexth + sum [ Li (Iai Ci + Pai) + X SbarTi taui ]
//
void cDynABNodeNOJ::_abBiasForce(cDynVector6& Pah, const cDynMatrix6& L, const cDynMatrix6& Ia, const cDynVector6& C, const cDynVector6& Pa)
{
    cDynVector6 tmpV61, tmpV62;

    tmpV61.multiply(Ia, C);
    tmpV61 += Pa;
    tmpV62.multiply(L, tmpV61);
    Pah += tmpV62;
}

//---------------------------------------------------------------------------
//
// Ii = Xc Ic Xtc
//    = [ RMRt, -RMRt rx; rx RMRt, RIRt - rx RMRt rx]
//    = [ M, -Mrx; rxM, RIRt - m rx rx]
//
void cDynABNodeNOJ::inertia(const double mass, const cDynVector3& centerOfMass, const cDynMatrix3& inertiaTensor)
{
    //  inertia = RIRt - mass()*rx*rx
    //  Ic = RIRt = inertia + mass()*rx*rx // change inertia
    cDynVector3 mr;
    mr.multiply(centerOfMass, mass);
    _I[0][1].cross(mr);  // mass()*rx
    _Ic.multiplyCross(_I[0][1], centerOfMass); // mass()*rx*rx
    _Ic += inertiaTensor; // Ic = RIRt

    _I[0][0].zero();
    _I[0][0][0][0] = _I[0][0][1][1] = _I[0][0][2][2] = mass;
    _I[0][1].negate(_I[0][1]);        // -mass()*rx
    _I[1][0].transpose(_I[0][1]);
    _I[1][1] = inertiaTensor;
}

//---------------------------------------------------------------------------
//
// FIX:  assuming Pi for cDynJoint = 0
//
// Pi = Xc (Wc X Ic Vc) - Ii (Wi X Vi)  
//    = Xc (Wc X [mi vc; Ic wc]) -Ii (Wi X Vi)
//
void cDynABNodeNOJ::biasForce(cDynVector6& P, const cDynVector6& V, const cDynVector3& WxV)
{
    // R = I_3;
    // Vc = Xtc Vi
    //    = Xt * V = [ Rt -Rtdx; 0 Rt ] [ v ; w ] = [ v-dxw ; w ]
    cDynVector6 tmpV6;
    double mass = _I[0][0][0][0];

    // mi vc = m ( v - dx w) = m v - m dx w = m v - I[1][0] w 
    // Ic wc = Ic w
    tmpV6[0].multiply(V[0], mass);
    tmpV6[1].multiply(_I[1][0], V[1]);

    tmpV6[0] -= tmpV6[1];			// mi vc
    tmpV6[1].multiply(_Ic, V[1]);	// Ic wc

    // Wc X tmpV61 = [ 0; w] X tmpV61 = [wx, 0 ; 0, wx] tmpV61
    P[0].crossMultiply(V[1], tmpV6[0]);
    P[1].crossMultiply(V[1], tmpV6[1]);

    // P = Xc (Wc X Ic Vc) = [ 1 , 0; rx , 1] [ p0 ; p1] = [p0; rx p0 + p1]
    tmpV6[0].multiply(_I[1][0], P[0]);

//	assert(mass > 0);
    
//	if (mass == 0)
//		mass = 1;

    tmpV6[0] *= 1/mass;
    P[1] += tmpV6[0];

    //  I (W X V) = [I00,I01;I10,I11]*[WxV;0] = [I00*WxV;I10*WxV]
    tmpV6[0].multiply(_I[0][0], WxV);
    tmpV6[1].multiply(_I[1][0], WxV);
    P -= tmpV6;
}

//---------------------------------------------------------------------------

double cDynABNodeNOJ::kineticEnergy(cDynVector6& V, const cDynVector6& Vh)
{
    cDynVector3 WxV, WhxVh;
    cDynVector6 IV;

    WhxVh.zero();

    velocity(V, WxV, Vh, WhxVh);
    
    IV.multiply(*I(), V);

    return V.dot(IV) / 2;
}

//---------------------------------------------------------------------------

double cDynABNodeNOJ::potentialEnergy(cDynVector3& g, const cDynVector3& gh, const cDynFrame& globalFrame, const double mass, const cDynVector3& centerOfMass)
{
    cDynVector3 h;
    cDynVector6 G;

    gravityForce(G, g, gh);

    h.multiply(globalFrame.rotation(), g);
    g = h;

    h.multiply(globalFrame.rotation(), centerOfMass);
    h += globalFrame.translation();
            
    return -mass * g.dot(h);	
}

//---------------------------------------------------------------------------

void cDynABNodeNOJ1::abImpulseInit()
{
    _joint->zero_Tau();

    flag(true);
}

//---------------------------------------------------------------------------

void cDynABNodeNOJn::abImpulseInit()
{
    for (int i = 0; i < noj(); i++)
        _joint[i]->zero_Tau();

    flag(true);
}

//---------------------------------------------------------------------------

void cDynABNodeNOJ1::updateLocalX(const cDynFrame& homeFrame)
{
    cDynTransform homeX;
    homeX.set(homeFrame);
    _joint->update_localX(homeX);
}

//---------------------------------------------------------------------------

void cDynABNodeNOJn::updateLocalX(const cDynFrame& homeFrame)
{
    cDynTransform homeX;
    homeX.set(homeFrame);
    _joint[0]->update_localX(homeX);

    homeX.identity();
    for (int i = 1; i < noj(); i++)
        _joint[i]->update_localX(homeX);

}

//---------------------------------------------------------------------------
//
// 0Jn = Jn = iXn^T Si = 0Xi^(-T) Si
// where 0Xi^(-T) = [ R rxR; 0 R ]
// since 0Rn = identity
// and  iXn^T = [(inv 0Xi) 0Xn]^T = (0Xn)^T * (0Xi)^(-T) = 0Xi ^(-T)
// since  0Xn = identity matrix    <--  {0} = {e}
//
void cDynABNodeNOJ1::globalJacobian(const cDynFrame& globalFrame)
{
    cDynTransform Xg;

    Xg.set(globalFrame);

    _joint->compute_Jg(Xg);
}

//---------------------------------------------------------------------------
//
// 0Xi = 0Xh hXi
// 0Xh = 0Xi hXi^-1
//
void cDynABNodeNOJn::globalJacobian(const cDynFrame& globalFrame)
{
    int i;
    cDynTransform Xg, Xgh;

    Xg.set(globalFrame);

    for (i = _noj - 1; i > 0; i--)
    {
        _joint[i]->compute_Jg(Xg);

        Xgh.multiplyInversed(Xg, _joint[i]->localX());
        Xg = Xgh;
    }
    _joint[i]->compute_Jg(Xg);
}

//---------------------------------------------------------------------------
//
// Gi = Xc Gc
// Gc = [ mi Rtci gi ; 0]
// gi = Rt gh
//
void cDynABNodeNOJ1::gravityForce(cDynVector6& G, cDynVector3& g, const cDynVector3& gh)
{
    g.transposedMultiply((_joint->localX().rotation()), gh);

    G[0].multiply(g, (*I())[0][0][0][0]);
    // [1,0;rx,1][g;0] = [g;rxg]
    G[1].multiply((*I())[1][0], g);
}

//---------------------------------------------------------------------------

void cDynABNodeNOJn::gravityForce(cDynVector6& G, cDynVector3& g, const cDynVector3& gh)
{
    cDynVector3 gp = gh;

    for (int i = 0; i < _noj; i++)
    {
        g.transposedMultiply((_joint[i]->localX().rotation()), gp);
        gp = g;
    }


    G[0].multiply(g, (*I())[0][0][0][0]);
    // [1,0;rx,1][g;0] = [g;rxg]
    G[1].multiply((*I())[1][0], g);
}

//---------------------------------------------------------------------------
//
// Vi = hXi^T Vh + Si dqi;
// xform = [R 0; dxR R]
// xformT = [ Rt -Rtdx; 0 Rt ]

// Ci = Wi X Vi - Xt (Wh X Vh) + Vi X Si dqi
// V X = [v0 ; v1] X = [ v1x , v0x ; 0 , v1x]
// WxV = [ 0 ; v1 ] x [ v0 ; v1 ]
//     = [ v1x , 0 ; 0 , v1x ] [ v0 ; v1 ] = [v1 x v0 ;v1 x v1] = [ v1 x v0 ; 0 ]
//     = [ wxv ; 0 ]
// Xt * WxV = [ Rt -Rtdx; 0 Rt ] [ wxv ; 0 ] = [ Rt WxV ; 0 ]
//
void cDynABNodeNOJ1::velocity(cDynVector6& V, cDynVector3& WxV, const cDynVector6& Vh, const cDynVector3& WhxVh)
{
    V.xformT(_joint->localX(), Vh);
    _joint->plusEq_SdQ(V);

    WxV.crossMultiply(V[1], V[0]);
    _joint->C()[0].transposedMultiply((_joint->localX().rotation()), WhxVh);
    _joint->C()[0].subtract(WxV, _joint->C()[0]);
    _joint->C()[1].zero();
    _joint->plusEq_V_X_SdQ(_joint->C(), V);
}

//---------------------------------------------------------------------------

void cDynABNodeNOJn::velocity(cDynVector6& V, cDynVector3& WxV, const cDynVector6& Vh, const cDynVector3& WhxVh)
{
    cDynVector6 Vp = Vh;
    cDynVector3 WpxVp = WhxVh;

    for (int i = 0; i < _noj; i++)
    {
        V.xformT(_joint[i]->localX(), Vp);
        _joint[i]->plusEq_SdQ(V);

        WxV.crossMultiply(V[1], V[0]);

        _joint[i]->C()[0].transposedMultiply((_joint[i]->localX().rotation()), WhxVh);
        _joint[i]->C()[0].subtract(WxV, _joint[i]->C()[0]);
        _joint[i]->C()[1].zero();
        _joint[i]->plusEq_V_X_SdQ(_joint[i]->C(), V);
    
        Vp = V;
        WpxVp = WxV;
    }
}

//---------------------------------------------------------------------------
//
// Dinv = inv(St Ia S)
// SbarT = Ia S Dinv
// hLi = hXi [ 1 - Si Sbari ]^T = hXi [1 - Sbari^T Si^T] = X - X SbarT St
// Lt = [1 - S Sbar] Xt
//
void cDynABNodeNOJ1::abInertiaDepend(cDynMatrix6& Iah, cDynVector6& Pah, cDynMatrix6& Ia, int propagate)
{
    _joint->minusEq_SdQ_damping(_joint->Pa(), Ia);

    _joint->compute_Dinv_and_SbarT(Ia);

    if (propagate)
    {
        _joint->L().set(_joint->localX());

        _joint->minusEq_X_SbarT_St(_joint->L(), _joint->localX());

        _abInertia(Iah, _joint->L(), Ia);
        _abBiasForce(Pah, _joint->L(), Ia, _joint->C(), _joint->Pa());
        _joint->plusEq_X_SbarT_Tau(Pah, _joint->localX());
    }
}

//---------------------------------------------------------------------------

void cDynABNodeNOJn::abInertiaDepend(cDynMatrix6& Iah, cDynVector6& Pah, cDynMatrix6& Ia, int propagate)
{
    cDynMatrix6 Iap;

    for (int i = _noj - 1; i > 0; i--)
    {
        _joint[i]->minusEq_SdQ_damping(_joint[i]->Pa(), Ia);
        _joint[i]->compute_Dinv_and_SbarT(Ia);

        _joint[i]->L().set(_joint[i]->localX());
        _joint[i]->minusEq_X_SbarT_St(_joint[i]->L(), _joint[i]->localX());

        _joint[i - 1]->Pa().zero();
        _abBiasForce(_joint[i - 1]->Pa(), _joint[i]->L(), Ia, _joint[i]->C(), _joint[i]->Pa());
        _joint[i]->plusEq_X_SbarT_Tau(_joint[i - 1]->Pa(), _joint[i]->localX());

        Iap.zero();
        _abInertia(Iap, _joint[i]->L(), Ia);

        Ia = Iap;
    }

    _joint[0]->minusEq_SdQ_damping(_joint[0]->Pa(), Ia);

    _joint[0]->compute_Dinv_and_SbarT(Ia);

    if (propagate)
    {
        _joint[0]->L().set(_joint[0]->localX());

        _joint[0]->minusEq_X_SbarT_St(_joint[0]->L(), _joint[0]->localX());

        _abInertia(Iah, _joint[0]->L(), Ia);
        _abBiasForce(Pah, _joint[0]->L(), Ia, _joint[0]->C(), _joint[0]->Pa());
        _joint[0]->plusEq_X_SbarT_Tau(Pah, _joint[0]->localX());
    }
}

//---------------------------------------------------------------------------

void cDynABNodeNOJ1::abInertiaDependConfig(cDynMatrix6& Iah, cDynMatrix6& Ia, int propagate)
{
    _joint->compute_Dinv_and_SbarT(Ia);

    if (propagate)
    {	
        _joint->L().set(_joint->localX());

        _joint->minusEq_X_SbarT_St(_joint->L(), _joint->localX());

        _abInertia(Iah, _joint->L(), Ia);
    }
}

//---------------------------------------------------------------------------

void cDynABNodeNOJn::abInertiaDependConfig(cDynMatrix6& Iah, cDynMatrix6& Ia, int propagate)
{
    cDynMatrix6 Iap;

    for (int i = _noj - 1; i > 0; i--)
    {
        _joint[i]->compute_Dinv_and_SbarT(Ia);

        _joint[i]->L().set(_joint[i]->localX());
        _joint[i]->minusEq_X_SbarT_St(_joint[i]->L(), _joint[i]->localX());

        Iap.zero();
        _abInertia(Iap, _joint[i]->L(), Ia);

        Ia = Iap;
    }

    _joint[0]->compute_Dinv_and_SbarT(Ia);

    _joint[0]->L().set(_joint[0]->localX());

    if (propagate)
    {
        _joint[0]->minusEq_X_SbarT_St(_joint[0]->L(), _joint[0]->localX());

        _abInertia(Iah, _joint[0]->L(), Ia);
    }
}

//---------------------------------------------------------------------------

void cDynABNodeNOJ1::abImpulse(cDynVector6& Yah, int propagate)
{
    if (propagate)
        Yah.multiply(_joint->L(), _joint->Pa());
}

//---------------------------------------------------------------------------

void cDynABNodeNOJn::abImpulse(cDynVector6& Yah, int propagate)
{
    for (int i = _noj - 1; i > 0; i--)
        _joint[i - 1]->Pa().multiply(_joint[i]->L(), _joint[i]->Pa());

    if (propagate)
        Yah.multiply(_joint[0]->L(), _joint[0]->Pa());
}

//---------------------------------------------------------------------------
//
// C = Ph = 0
// Pah = - Fexth + sum [ Li (Pai) + X SbarTi taui ]
//
void cDynABNodeNOJ1::abBiasForceConfig(cDynVector6& Pah, int propagate)
{
    // store fwdIn() in dq_ to be used by JointAcceleration(true)
    //_joint[0]->read_dQ(fwdIn, fwdInSphere);

    if (propagate)
    {
        Pah.multiply(_joint->L(), _joint->Pa());
        _joint->plusEq_X_SbarT_Tau(Pah, _joint->localX());
    }
}

//---------------------------------------------------------------------------

void cDynABNodeNOJn::abBiasForceConfig(cDynVector6& Pah, int propagate)
{
    // store fwdIn() in dq_ to be used by JointAcceleration(true)
    //_joint[0]->read_dQ(fwdIn, fwdInSphere);

    for (int i = _noj - 1; i > 0; i--)
    {
        _joint[i - 1]->Pa().multiply(_joint[i]->L(), _joint[i]->Pa());
        _joint[i]->plusEq_X_SbarT_Tau(_joint[i - 1]->Pa(), _joint[i]->localX());
    }
    if (propagate)
    {
        Pah.multiply(_joint[0]->L(), _joint[0]->Pa());
        _joint[0]->plusEq_X_SbarT_Tau(Pah, _joint[0]->localX());
    }
}

//---------------------------------------------------------------------------
//
// ddQ = Dinv*(tau - St*Pa) - Sbar*(X Ah + Ci)
// Ai = (hXi^T Ah + Ci) + Si ddqi;
//
void cDynABNodeNOJ1::acceleration(cDynVector6& A, const cDynVector6& Ah)
{
    A.xformT(_joint->localX(), Ah);
    A += _joint->C();

    _joint->compute_ddQ(_joint->Pa(), A);
    _joint->plusEq_SddQ(A);
}

//---------------------------------------------------------------------------

void cDynABNodeNOJn::acceleration(cDynVector6& A, const cDynVector6& Ah)
{
    cDynVector6 Ap = Ah;

    for (int i = 0; i < _noj; i++)
    {
        A.xformT(_joint[i]->localX(), Ap);
        A += _joint[i]->C();
        _joint[i]->compute_ddQ(_joint[i]->Pa(), A);
        _joint[i]->plusEq_SddQ(A);
        Ap = A;
    }
}

//---------------------------------------------------------------------------

void cDynABNodeNOJ1::accelerationOnly(cDynVector6& A, const cDynVector6& Ah)
{
    A.xformT(_joint->localX(), Ah);
    A += _joint->C();
    _joint->plusEq_SddQ(A);
}

//---------------------------------------------------------------------------

void cDynABNodeNOJn::accelerationOnly(cDynVector6& A, const cDynVector6& Ah)
{
    cDynVector6 Ap = Ah;

    for (int i = 0; i < _noj; i++)
    {
        A.xformT(_joint[i]->localX(), Ap);
        A += _joint[i]->C();
        _joint[i]->plusEq_SddQ(A);
        Ap = A;
    }
}

//---------------------------------------------------------------------------

void cDynABNodeNOJ1::force(cDynVector6& Fh, int propagate)
{
    // no damping
    if (propagate)
    {
        cDynVector6 tmpV;
        tmpV.xform(_joint->localX(), _joint->Pa());
        Fh += tmpV;
    }
    _joint->compute_Tau(_joint->Pa());
}

//---------------------------------------------------------------------------

void cDynABNodeNOJn::force(cDynVector6& Fh, int propagate)
{
    for (int i = _noj - 1; i > 0; i--)
    {
        _joint[i - 1]->Pa().xform(_joint[i]->localX(), _joint[i]->Pa());
        _joint[i]->compute_Tau(_joint[i]->Pa());
    }

    if (propagate)
    {
        cDynVector6 tmpV;
        tmpV.xform(_joint[0]->localX(), _joint[0]->Pa());
        Fh += tmpV;
    }
    _joint[0]->compute_Tau(_joint[0]->Pa());
}

//---------------------------------------------------------------------------
//
// d dQ = Dinv*(- St*Pa) - Sbar*(X dVh)
// dVi = (hXi^T dVh) + Si d dqi;
//
// ddQ = Dinv*(tau - St*Pa) - Sbar*(X Ah)
// Ai = (hXi^T Ah) + Si ddqi;
//
void cDynABNodeNOJ1::velocityDelta(cDynVector6& dV, const cDynVector6& dVh)
{
    dV.xformT(_joint->localX(), dVh);
    if (flag())
    {
        _joint->compute_ddQ(_joint->Pa(), dV);
    }
    else
    {
        _joint->compute_ddQ(dV);
    }
    _joint->plusEq_SddQ(dV);

    flag(0);
}

//---------------------------------------------------------------------------

void cDynABNodeNOJn::velocityDelta(cDynVector6& dV, const cDynVector6& dVh)
{
    cDynVector6 dVp = dVh;

    for (int i = 0; i < _noj; i++)
    {
        dV.xformT(_joint[i]->localX(), dVp);
        if (flag())
        {
            _joint[i]->compute_ddQ(_joint[i]->Pa(), dV);		
        }
        else
        {
            _joint[i]->compute_ddQ(dV);
        }
        _joint[i]->plusEq_SddQ(dV);
        dVp = dV;
    }

    flag(0);
}

//---------------------------------------------------------------------------



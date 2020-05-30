//===========================================================================
/*
    This file is part of the dynamics library.
    Copyright (C) 2014, Artificial Intelligence Laboratory,
    Stanford University. All rights reserved.

    \version   $MAJOR.$MINOR.$RELEASE $Rev: 1242 $
*/
//===========================================================================

//---------------------------------------------------------------------------
#include "dynamics/CDynABJoint.h"
//---------------------------------------------------------------------------
#include <assert.h>
//---------------------------------------------------------------------------

void cDynABJointDOF1::update_localX(const cDynTransform& home)
{
    cDynTransform local;

    if (_S[1].dot(_S[1]) > CDYN_QUATERNION_EPSILON)
    {
        cDynQuaternion tmpQ;
        tmpQ.set(_S[1], _Q);
        tmpQ.normalize();
        local.set(tmpQ);	
    }
    else
        local.identity();

    cDynVector3 tmpV3;
    tmpV3.multiply(_S[0], _Q);
    local.set(tmpV3);
    
    localX().multiply(home, local);
}


//---------------------------------------------------------------------------
//
// Vi = hXi^T Vh + Si dqi;
// xform = [R 0; dxR R]
// xformT = [ Rt -Rtdx; 0 Rt ]
//
void cDynABJointDOF1::plusEq_SdQ(cDynVector6& V)
{
    cDynVector6 tmpV6;
    tmpV6.multiply(_S, _dQ);
    V += tmpV6;
}


//---------------------------------------------------------------------------
//
// Ci = Wi X Vi - Xt (Wh X Vh) + Vi X Si dqi
// V X = [v0 ; v1] X = [ v1x , v0x ; 0 , v1x]
// WxV = [ 0 ; v1 ] x [ v0 ; v1 ]
//     = [ v1x , 0 ; 0 , v1x ] [ v0 ; v1 ] = [v1 x v0 ;v1 x v1] = [ v1 x v0 ; 0 ]
//     = [ wxv ; 0 ]
// Xt * WxV = [ Rt -Rtdx; 0 Rt ] [ wxv ; 0 ] = [ Rt WxV ; 0 ]
//
void cDynABJointDOF1::plusEq_V_X_SdQ(cDynVector6& C, const cDynVector6& V)
{
    cDynVector6 tmpV6;
    tmpV6.crossMultiply(V, _S);
    tmpV6 *= _dQ;
    C += tmpV6;
}

//---------------------------------------------------------------------------
//
// Dinv = inv(St Ia S)
// SbarT = Ia S Dinv
// hLi = hXi [ 1 - Si Sbari ]^T = hXi [1 - Sbari^T Si^T] = X - X SbarT St
// Lt = [1 - S Sbar] Xt
//
void cDynABJointDOF1::compute_Dinv_and_SbarT(const cDynMatrix6& Ia)
{
    cDynVector6 IaS;
    IaS.multiply(Ia, _S);
    _Dinv = _S.dot(IaS) + _motorInertia;
    _Dinv = 1/_Dinv;

    _SbarT.multiply(IaS, _Dinv);
}

//---------------------------------------------------------------------------

void cDynABJointDOF1::minusEq_X_SbarT_St(cDynMatrix6& L, const cDynTransform& localX)
{
    cDynVector6 tmpV6;
    tmpV6.xform(localX, _SbarT);
    cDynMatrix6 tmpM6;
    tmpM6.multiplyTransposed(tmpV6, _S);
    L -= tmpM6;
}

//---------------------------------------------------------------------------
//
// Pah = Ph - Fexth + sum [ Li (Iai Ci + Pai) + X SbarTi taui ]
//
void cDynABJointDOF1::plusEq_X_SbarT_Tau(cDynVector6& Pah, const cDynTransform& localX)
{
    cDynVector6 tmpV6;
    tmpV6.xform(localX, _SbarT);
    tmpV6 *= _Tau;
    Pah += tmpV6;
}

//---------------------------------------------------------------------------
//
void cDynABJointDOF1::compute_Tau(const cDynVector6& F)
{
    _Tau = _S.dot(F) + _ddQ * _motorInertia;
}

//---------------------------------------------------------------------------
//
// ddQ = Dinv*(tau - St*Pa) - Sbar*(X Ah + Ci)
// Ai = (hXi^T Ah + Ci) + Si ddqi;
//
void cDynABJointDOF1::compute_ddQ(const cDynVector6& Pa, const cDynVector6& XAh_C)
{
    _ddQ = _Dinv * (_Tau - _S.dot(Pa)) - _SbarT.dot(XAh_C);
}

//---------------------------------------------------------------------------

void cDynABJointDOF1::compute_ddQ(const cDynVector6& XAh_C)
{
    _ddQ = -_SbarT.dot(XAh_C);
}

//---------------------------------------------------------------------------

void cDynABJointDOF1::plusEq_SddQ(cDynVector6& A)
{
    cDynVector6 tmpV6;
    tmpV6.multiply(_S, _ddQ);
    A += tmpV6;
}

//---------------------------------------------------------------------------

void cDynABJointDOF1::minusEq_SdQ_damping(cDynVector6& B, const cDynMatrix6& Ia)
{
    cDynVector6 tmpV;
    tmpV.multiply(Ia, _S);
    tmpV *= _dQ * (- _damping);
    B -= tmpV;
}

//---------------------------------------------------------------------------
//
// 0Jn = Jn = iXn^T Si = 0Xi^(-T) Si
// where 0Xi^(-T) = [ R rxR; 0 R ]
//
void cDynABJointDOF1::compute_Jg(const cDynTransform &globalX)
{ 
    _Jg.xformInvT(globalX, _S);
}

//---------------------------------------------------------------------------

void cDynABJointSpherical::update_localX(const cDynTransform& home)
{
    cDynTransform local;

    local.set(_Q);
    local.set(0, 0, 0);

    localX().multiply(home, local);
}

//---------------------------------------------------------------------------
//
// Vi = hXi^T Vh + Si dqi;
// xform = [R 0; dxR R]
// xformT = [ Rt -Rtdx; 0 Rt ]
//
void cDynABJointSpherical::plusEq_SdQ(cDynVector6& V)
{
    V[1] += _dQ;
}

//---------------------------------------------------------------------------
//
// Ci = Wi X Vi - Xt (Wh X Vh) + Vi X Si dqi
// V X = [v0 ; v1] X = [ v1x , v0x ; 0 , v1x]
// WxV = [ 0 ; v1 ] x [ v0 ; v1 ]
//     = [ v1x , 0 ; 0 , v1x ] [ v0 ; v1 ] = [v1 x v0 ;v1 x v1] = [ v1 x v0 ; 0 ]
//     = [ wxv ; 0 ]
// Xt * WxV = [ Rt -Rtdx; 0 Rt ] [ wxv ; 0 ] = [ Rt WxV ; 0 ]
// V X [0 ; dq] = [ v1x , v0x ; 0 , v1x] [0 ; dq] = [ v0 x dq ; v1 x dq]
//
void cDynABJointSpherical::plusEq_V_X_SdQ(cDynVector6& C, const cDynVector6& V)
{
    cDynVector3 tmpV3;
    tmpV3.crossMultiply(V[0], _dQ);
    C[0] += tmpV3;
    tmpV3.crossMultiply(V[1], _dQ);
    C[1] += tmpV3;
}

//---------------------------------------------------------------------------
//
// Dinv = inv(St Ia S)
// SbarT = Ia S Dinv
// hLi = hXi [ 1 - Si Sbari ]^T = hXi [1 - Sbari^T Si^T] = X - X SbarT St
// Lt = [1 - S Sbar] Xt
//
void cDynABJointSpherical::compute_Dinv_and_SbarT(const cDynMatrix6& Ia)
{
    cDynMatrix3 I = Ia[1][1];
    I[0][0] += _motorInertia;
    I[1][1] += _motorInertia;
    I[2][2] += _motorInertia;
    _Dinv.inverseDetSPD(I);

    _SbarT[0].multiply(Ia[0][1], _Dinv);
    _SbarT[1].multiply(Ia[1][1], _Dinv);
}

//---------------------------------------------------------------------------
//
// Xform = [R 0; dxR R]
// L_ = X[1 - SbarT St] = X - X SbarT St
// SbarT St = [0, SbarT0; 0 SbarT1]
// X SbarT St = [R ,0; dxR, R][0, SbarT0; 0 SbarT1]=[0,R SbarT0;0, dxR SbarT0 + R SbarT1]
// L = X - X SbarT St = [ R , -R SbarT0; dxR , R - (dxR SbarT0 + R SbarT1) ]
//                                            = R (1 - SbarT1) - dxR SbarT0
//
void cDynABJointSpherical::minusEq_X_SbarT_St(cDynMatrix6& L, const cDynTransform& localX)
{
    cDynMatrix3 tmpM31, tmpM32;
    tmpM31.multiply(localX.rotation(), _SbarT[0]);
    L[0][1] -= tmpM31;
    tmpM32.multiply(localX.rotation(), _SbarT[1]);
    L[1][1] -= tmpM32;
    tmpM32.crossMultiply(localX.translation(), tmpM31);
    L[1][1] -= tmpM32;
}

//---------------------------------------------------------------------------
//
// Pah = Ph - Fexth + sum [ Li (Iai Ci + Pai) + X SbarTi taui ]
//
void cDynABJointSpherical::plusEq_X_SbarT_Tau(cDynVector6& Pah, const cDynTransform& localX)
{
    cDynVector6 tmpV61, tmpV62;
    tmpV61[0].multiply(_SbarT[0], _Tau);
    tmpV61[1].multiply(_SbarT[1], _Tau);
    tmpV62.xform(localX, tmpV61);
    Pah += tmpV62;
}

//---------------------------------------------------------------------------

void cDynABJointSpherical::compute_Tau(const cDynVector6& F)
{
    _Tau.multiply(_ddQ, _motorInertia);
    _Tau += F[1];
}

//---------------------------------------------------------------------------
//
// ddQ = Dinv*(tau - St*Pa) - Sbar*(X Ah + Ci)
// Ai = (hXi^T Ah + Ci) + Si ddqi;
//
void cDynABJointSpherical::compute_ddQ(const cDynVector6& Pa, const cDynVector6& XAh_C)
{
    cDynVector3 tmpV3;
    // St Pa = Pa[1]
    tmpV3.subtract(_Tau, Pa[1]);
    _ddQ.multiply(_Dinv, tmpV3);  // ddq = Dinv * (tau - Pa)
    // Sbar*(X Ah + Ci) = (SbarT)^t * (X Ah + Ci)
    tmpV3.transposedMultiply(_SbarT[0], XAh_C[0]);
    _ddQ -= tmpV3;
    tmpV3.transposedMultiply(_SbarT[1], XAh_C[1]);
    _ddQ -= tmpV3;
}

//---------------------------------------------------------------------------
//
void cDynABJointSpherical::compute_ddQ(const cDynVector6& XAh_C)
{
    cDynVector3 tmpV3;

    // Sbar*(X Ah + Ci) = (SbarT)^t * (X Ah + Ci)
    tmpV3.transposedMultiply(_SbarT[0], XAh_C[0]);
    _ddQ.negate(tmpV3);
    tmpV3.transposedMultiply(_SbarT[1], XAh_C[1]);
    _ddQ -= tmpV3;
}

//---------------------------------------------------------------------------
//
void cDynABJointSpherical::plusEq_SddQ(cDynVector6& A)
{
    A[1] += _ddQ;
}

//---------------------------------------------------------------------------
//
void cDynABJointSpherical::minusEq_SdQ_damping(cDynVector6& B, const cDynMatrix6& Ia)
{
    cDynVector3 tmpV;
    tmpV.multiply(Ia[1][1], _dQ);
    tmpV *= (- _damping);
    B[1] -= tmpV;
}

//---------------------------------------------------------------------------
//
// 0Jn = Jn = iXn^T Si = 0Xi^(-T) Si
// where 0Xi^(-T) = [ R rxR; 0 R ]
//
void cDynABJointSpherical::compute_Jg(const cDynTransform &globalX) 
{ 
    _Jg[0].crossMultiply(globalX.translation(), globalX.rotation());
    _Jg[1] = globalX.rotation(); 
}

//---------------------------------------------------------------------------

void cDynABJointFree::update_localX(const cDynTransform& home)
{
    cDynTransform local;

    local.set(_Q);

    localX().multiply(home, local);
}

//---------------------------------------------------------------------------
//
// Vi = hXi^T Vh + Si dqi;
// xform = [R 0; dxR R]
// xformT = [ Rt -Rtdx; 0 Rt ]
//
void cDynABJointFree::plusEq_SdQ(cDynVector6& V)
{
    V += _dQ;
}

//---------------------------------------------------------------------------
//
// Ci = Wi X Vi - Xt (Wh X Vh) + Vi X Si dqi
// V X = [v0 ; v1] X = [ v1x , v0x ; 0 , v1x]
// WxV = [ 0 ; v1 ] x [ v0 ; v1 ]
//     = [ v1x , 0 ; 0 , v1x ] [ v0 ; v1 ] = [v1 x v0 ;v1 x v1] = [ v1 x v0 ; 0 ]
//     = [ wxv ; 0 ]
// Xt * WxV = [ Rt -Rtdx; 0 Rt ] [ wxv ; 0 ] = [ Rt WxV ; 0 ]
//
void cDynABJointFree::plusEq_V_X_SdQ(cDynVector6& C, const cDynVector6& V)
{
//	cDynVector6 VxSdQ;
//	VxSdQ.crossMultiply(V, _dQ);
//	C += VxSdQ;
    C.zero();
}

//---------------------------------------------------------------------------
//
// Dinv = inv(St Ia S)
// SbarT = Ia S Dinv
// hLi = hXi [ 1 - Si Sbari ]^T = hXi [1 - Sbari^T Si^T] = X - X SbarT St
// Lt = [1 - S Sbar] Xt
//
void cDynABJointFree::compute_Dinv_and_SbarT(const cDynMatrix6& Ia)
{
    cDynMatrix6 I = Ia;
    for (int i = 0; i < 6; i++)
        I.elementAt(i, i) += _motorInertia;
    // Dinv = inverse Ia
    _Dinv.luDecomp(I);
    //_Dinv.luDecompSPD(Ia);
    // SbarT = identity
}

//---------------------------------------------------------------------------

void cDynABJointFree::minusEq_X_SbarT_St(cDynMatrix6& L, const cDynTransform& localX)
{
    L.zero();
}

//---------------------------------------------------------------------------
//
// Pah = Ph - Fexth + sum [ Li (Iai Ci + Pai) + X SbarTi taui ]
//
void cDynABJointFree::plusEq_X_SbarT_Tau(cDynVector6& Pah, const cDynTransform& localX)
{
    assert(false);

//	cDynVector6 tmpV;
//	tmpV.xform(localX, _Tau);
//	Pah += tmpV;
}

//---------------------------------------------------------------------------

void cDynABJointFree::compute_Tau(const cDynVector6& F)
{
    _Tau.multiply(_ddQ, _motorInertia);
    _Tau += F;
}

//---------------------------------------------------------------------------
//
// ddQ = Dinv*(tau - St*Pa) - Sbar*(X Ah + Ci)
// Ai = (hXi^T Ah + Ci) + Si ddqi;
//
void cDynABJointFree::compute_ddQ(const cDynVector6& Pa, const cDynVector6& XAh_C)
{
    // Dinv * (tau - Pa) = ddQ + XAh + C
    // (tau - Pa) = Ia (ddQ + XAh + C)
    //  let x = ddQ + XAh + C
    // solve for x and then ddQ = x - (XAh + C)
    cDynVector6 tmpV;
    tmpV.subtract(_Tau, Pa);
    _ddQ.backSub(_Dinv, tmpV);
    //_ddQ.backSubSPD(_Dinv, tmpV);
    // note that XAh + C = 0
//	_ddQ -= XAh_C;
}

//---------------------------------------------------------------------------

void cDynABJointFree::compute_ddQ(const cDynVector6& XAh_C)
{
    assert(false);

//	_ddQ.negate(XAh_C);
}

//---------------------------------------------------------------------------

void cDynABJointFree::plusEq_SddQ(cDynVector6& A)
{
    A += _ddQ;
}

//---------------------------------------------------------------------------

void cDynABJointFree::minusEq_SdQ_damping(cDynVector6& B, const cDynMatrix6& Ia)
{
    cDynVector6 tmpV6;
    tmpV6.multiply(Ia, _dQ);
    tmpV6 *= (- _damping);
    B -= tmpV6; 
}

//---------------------------------------------------------------------------
//
// 0Jn = Jn = iXn^T Si = 0Xi^(-T) Si
// where 0Xi^(-T) = [ R rxR; 0 R ]
//
void cDynABJointFree::compute_Jg(const cDynTransform &globalX)
{ 
    _Jg[0][0] = globalX.rotation(); 
    _Jg[0][1].crossMultiply(globalX.translation(), globalX.rotation());
    _Jg[1][0].zero();
    _Jg[1][1] = globalX.rotation(); 
}

//---------------------------------------------------------------------------



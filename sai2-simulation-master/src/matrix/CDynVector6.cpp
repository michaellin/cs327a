//===========================================================================
/*
    This file is part of the dynamics library.
    Copyright (C) 2014, Artificial Intelligence Laboratory,
    Stanford University. All rights reserved.

    \version   $MAJOR.$MINOR.$RELEASE $Rev: 1242 $
*/
//===========================================================================

//---------------------------------------------------------------------------
#include "matrix/CDynVector6.h"
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
/* 
    y = LU x
    find x given LU and y 
    diag(U) = [1 1 1 ... 1]
 */
//---------------------------------------------------------------------------
void cDynVector6::backSub(const cDynMatrix6 &lu, const cDynVector6 &y)
{
    //
    // LU x = y -> x 
    // L d  = y -> fw sub  
    // U x = d -> bw sub  
    //

    int i, j;
    double d[6];

    /// forward
    for (i = 0; i < 6; i++)
    {
        d[i] = y.elementAt(i);
        for (j = 0; j < i; j++)
            d[i] -= lu.elementAt(i, j) * d[j];
        d[i] /= lu.elementAt(i, i);
    }
    
    // backward
    for (i = 5; i >= 0; i--)
    {
        elementAt(i) = d[i];
        for (j = i + 1; j < 6; j++)
            elementAt(i) -= lu.elementAt(i, j) * elementAt(j);
    }
}

//---------------------------------------------------------------------------
/*
    y = LU x
    y = LL' x
    find x given LU and y 
    (symmetric and positive definite)
 */
//---------------------------------------------------------------------------
void cDynVector6::backSubSPD(const cDynMatrix6 &lu, const cDynVector6 &y)
{
    //
    // U = L'  
    // LU = LL'  
    // LU x = y -> x  
    // L d  = y -> fw sub  
    // U x = d -> bw sub  
    //

    int i, j;
    double d[6];

    // forward
    for (i = 0; i < 6; i++)
    {
        d[i] = y.elementAt(i);
        for (j = 0; j < i; j++)
            d[i] -= lu.elementAt(i, j) * d[j];
        d[i] /= lu.elementAt(i, i);
    }
    // backward
    for (i = 5; i >= 0; i--)
    {
        elementAt(i) = d[i];
        for (j = i + 1; j < 6; j++)
            elementAt(i) -= lu.elementAt(j, i) * elementAt(j);
        elementAt(i) /= lu.elementAt(i, i);
    }
}

//---------------------------------------------------------------------------

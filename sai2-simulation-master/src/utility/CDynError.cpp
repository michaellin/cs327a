//===========================================================================
/*
    This file is part of the dynamics library.
    Copyright (C) 2014, Artificial Intelligence Laboratory,
    Stanford University. All rights reserved.

    \version   $MAJOR.$MINOR.$RELEASE $Rev: 1242 $
*/
//===========================================================================

//---------------------------------------------------------------------------
#include <stdlib.h>
#include "utility/CDynLogger.h"
#include "utility/CDynError.h"
//---------------------------------------------------------------------------
cbDeError g_cbDynError = NULL;
//---------------------------------------------------------------------------

void cDynErrorSetCallback(cbDeError f)
{
    g_cbDynError = f;
}

//---------------------------------------------------------------------------

void cDynError(cDynErrorType error, void* data)
{
    if (g_cbDynError)
    {
        g_cbDynError(error, data);
    }
    else
    {
        cDynPrintf("ERROR (%d)\n",(int)error);
    }
}

//---------------------------------------------------------------------------
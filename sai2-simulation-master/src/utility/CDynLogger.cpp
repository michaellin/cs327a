//===========================================================================
/*
    This file is part of the dynamics library.
    Copyright (C) 2014, Artificial Intelligence Laboratory,
    Stanford University. All rights reserved.

    \version   $MAJOR.$MINOR.$RELEASE $Rev: 1242 $
*/
//===========================================================================

//---------------------------------------------------------------------------
#include <stdio.h> // NULL
#include <stdarg.h>
//---------------------------------------------------------------------------
#if defined(WIN32) | defined(WIN64)
#ifndef BORLAND
#include <windows.h> // OutputDebugString
#endif
#endif
//---------------------------------------------------------------------------
#include "utility/CDynLogger.h"
//---------------------------------------------------------------------------
#if defined(WIN32) | defined(WIN64)
#ifndef BORLAND
#define vsnprintf _vsnprintf
#endif
#endif
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
// DEFINES
//---------------------------------------------------------------------------

#define	PRINTF_BUFSIZE	4096

//---------------------------------------------------------------------------
// STATICS
//---------------------------------------------------------------------------

cDynLogger* cDynLogger::logger = NULL;

//---------------------------------------------------------------------------

void cDynLoggerOutputStd::Log(char *msg)
{
    printf("%s", msg);
    fflush(stdout);
}

//---------------------------------------------------------------------------

#if defined(WIN32) | defined(WIN64)
#ifndef BORLAND
void cDynLoggerOutputWinDebug::Log(char *msg)
{
    OutputDebugString(msg);
}
#endif
#endif

//---------------------------------------------------------------------------

cDynLoggerOutputFile::cDynLoggerOutputFile(char* filename)
{
    logfile = fopen(filename, "w");
}

//---------------------------------------------------------------------------

cDynLoggerOutputFile::~cDynLoggerOutputFile()
{
    if (logfile) fclose(logfile);
}

//---------------------------------------------------------------------------

void cDynLoggerOutputFile::Log(char *msg)
{
    if (!logfile) return;
    fprintf(logfile, "%s", msg);
}

//---------------------------------------------------------------------------

cDynLogger::cDynLogger()
{
    output = NULL;
}

//---------------------------------------------------------------------------

cDynLogger::~cDynLogger()
{
    if (output) delete output;
}

//---------------------------------------------------------------------------

void cDynLogger::Initialize()
{
    if (!logger) 
    {
        logger = new cDynLogger;
    }
}

//---------------------------------------------------------------------------

void cDynLogger::Shutdown()
{
    if (logger)
        delete logger;
}

//---------------------------------------------------------------------------

void cDynLogger::Format(char *buf, int type, char* category, char *msg)
{
    char *prefix = NULL;
    char str[12];

    if (type == cDynLogger::_DFL_ERROR) prefix = "Error";
    else if (type == cDynLogger::_DFL_WARN)  prefix = "Warning";
    else if (type == cDynLogger::_DFL_INFO)  prefix = "Info";
    else if (type == cDynLogger::_DFL_DEBUG) prefix = "Debug";
    else 
    {
        sprintf(str, "%d", type);
        prefix = str;
    }

    if (category)
        sprintf(buf, "%s(%s): %s\n", prefix, category, msg);
    else
        sprintf(buf, "%s: %s\n", prefix, msg);
}

//---------------------------------------------------------------------------

void cDynLogger::Log(int type, char* category, char* format, ...)
{
    va_list argptr;
    char buf[PRINTF_BUFSIZE], msg[PRINTF_BUFSIZE];
    cDynLoggerOutput *pout;

    if (!logger || !logger->output) return;

    va_start(argptr, format);
#ifdef PS2
    vsprintf(msg, format, argptr);
#else
    vsnprintf(msg, PRINTF_BUFSIZE, format, argptr);
#endif
    va_end(argptr);

    Format(buf, type, category, msg);

    pout = logger->output;
    while (pout)
    {
        pout->Log(buf);
        pout = pout->next;
    }
}

//---------------------------------------------------------------------------

void cDynLogger::Printf(char* format, ...)
{
    va_list argptr;
    char buf[PRINTF_BUFSIZE];
    cDynLoggerOutput *pout;

    if (!logger || !logger->output) return;

    va_start(argptr, format);
#ifdef PS2
    vsprintf(buf, format, argptr);
#else
    vsnprintf(buf, PRINTF_BUFSIZE, format, argptr);
#endif
    va_end(argptr);

    pout = logger->output;
    while (pout)
    {
        pout->Log(buf);
        pout = pout->next;
    }
}

//---------------------------------------------------------------------------

void cDynLogger::AddOutput(cDynLoggerOutput *out)
{
    Initialize();

    out->next = logger->output;
    logger->output = out;
}

//---------------------------------------------------------------------------

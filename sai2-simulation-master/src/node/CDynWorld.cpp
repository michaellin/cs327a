//===========================================================================
/*
    This file is part of the dynamics library.
    Copyright (C) 2014, Artificial Intelligence Laboratory,
    Stanford University. All rights reserved.

    \version   $MAJOR.$MINOR.$RELEASE $Rev: 1242 $
*/
//===========================================================================

//---------------------------------------------------------------------------
#include "node/CDynWorld.h"
//---------------------------------------------------------------------------
//#include "utility/CDynTimeInterval.h"
#include "utility/CDynLogger.h"
#include "utility/CDynError.h"
#include "object/CDynObject.h"
#include "dynamics/CDynDynamics.h"
#include "node/CDynContact.h"
#include "node/CDynContactPointCache.h"
#include "node/CDynPrimPair.h"
#include "node/CDynConstraintList.h"
//---------------------------------------------------------------------------
// global variable used for debugging
long g_worldIter=0;
static char* (*cDynObjectNameFunction)(cDynObject* world, void* arg)=NULL;
//---------------------------------------------------------------------------
cDynVector3 cVector3Zero;
cDynMatrix3 cMatrix3Zero;
//---------------------------------------------------------------------------
// bool g_noBackup=false;
//---------------------------------------------------------------------------

void cDynIntegrateCallback(cDynTimeEvent* event, void *arg)
{
    cDynBaseNode* base=(cDynBaseNode *)arg;
    cDynTime inc=event->inc();

    //cDynPrintf("INTEGRATE %5.9f\n",event->time());
    base->integrate(event->time(),inc);
}

//---------------------------------------------------------------------------

void cDynObjectNameCallback(char* (*f)(cDynObject* obj, void* arg))
{
    cDynObjectNameFunction=f;
}

//---------------------------------------------------------------------------

char *cDynObjectName(cDynObject* obj)
{
    if (cDynObjectNameFunction)
        return(cDynObjectNameFunction(obj,obj->data()));
    else return("none");
}

//---------------------------------------------------------------------------

static void cDynCallbackCollisionTest(cDynTimeEvent* event, void* arg)
{
    cDynWorld* world=(cDynWorld *)arg;
    
    //cDynPrintf("COLLEVENT %5.9f\n",event->time());
    world->collisionTest(event->time());
}


//---------------------------------------------------------------------------

cDynWorld::cDynWorld()
{
    time_=0.0;

    backupTime_=0.0;
    backupLimit_=0.0;

    firstBase_=lastBase_=NULL;

    collisionEvent_=new cDynTimeEvent(time_, -1, 1.0/30.0, CDYN_TIMEEVENT_FOREVER,&cDynCallbackCollisionTest,(void *)this);
    timer_.insert(collisionEvent_);
    
    // gravity.set(0.0,0.0,-9.81);
    gravity.set(0.0,0.0,-10.0);
    
    callback_=NULL;

    cVector3Zero.zero();
    cMatrix3Zero.zero();
}

//---------------------------------------------------------------------------

cDynBaseNode* cDynWorld::insert(cDynObject* obj, char* data)
{
    cDynBaseNode *baseNode=new cDynBaseNode(this);
    baseNode->link(obj);
    baseNode->data(data);
    baseNode->insert(&firstBase_,&lastBase_);
    return(baseNode);
}

//---------------------------------------------------------------------------

void cDynWorld::remove(cDynBaseNode* base)
{
    assert(base->world() == this);

    base->remove(&firstBase_,&lastBase_);
}

//---------------------------------------------------------------------------

void cDynWorld::backup(const cDynTime &time)
{
    cDynTime t=time;

    // negative means relative time
    if (t < 0.0) t += time_;	

    // time is later than current time
    if (t > timer_.top()->time()) return;

    // backup is beyond limits
    if (t < backupTime_) t=backupTime_; 

    for (cDynBaseNode* b=firstBase_;b != NULL; b=b->next_) 
    {
        b->backup(time);
    }
    time_=t;
}

//---------------------------------------------------------------------------

void cDynWorld::collisionTest(const cDynTime &tc)
{
    static int flag = 0;
    g_worldIter++;
    //time_=tc;
    if (time_ == tc) return;

    if (cDynDebugLevel >= 1)
        cDynPrintf("TEST %5.5f to %5.5f\n",time_,tc);

    for (cDynBaseNode* b=firstBase_;b != NULL; b=b->next_) 
    {
        if (b->status() == CDYN_ACTIVE) 
        {
            cDynTime time=tc;

            b->contact()->updatePersistent(time);
            b->checkInvalids(time);
            time=b->checkJointLimits(time);
            b->root()->collisionUpdateTree(time_,time);
        }
    }

    // check user defined constraints
    constraints_.check(time_,tc);
    
    // check for collision between objects
    if (!collision_.checkContacts()) 
    { 
        // an invalid constraint was found reschedule collision test
        collisionEvent_->time(tc);
        return;
    }
#if 0
    if (flag) 
    {
        cDynPrintf("Before\n");
        for (cDynBaseNode* b=firstBase_;b != NULL; b=b->next_) 
        {
            b->rootNode()->displaySmallTree(NULL);
        }
    }
#endif
    //backup everyone to their contact time
    for (cDynBaseNode* b=firstBase_;b != NULL; b=b->next_) 
    {
        if (b->contact()->n())
            b->backup(b->contact()->time());
    }
    // for each group update
    for (cDynBaseNode* b=firstBase_;b != NULL; b=b->next_) 
    {
        if (!b->contact()->setup() && b->contact()->n()) 
        {
#ifdef C_DYN_DEBUG_OPS
    cDynDebugEnergyReset();
#endif
            g_dynTheContactPointCache.reset();
            int m=b->contact()->setup(0);

            //rid += m;
            if (!b->contact()->resolve(m)) 
            {
                // an action has changed the state of the world
                // most likely a CDYN_SLEEP base has been made 
                // ACTIVE. reschedule collisionTest 
                collisionEvent_->time(b->contact()->time());
                return;
            }

#ifdef C_DYN_DEBUG_OPS
            double deltaKE=cDynDebugTotal(1)-cDynDebugTotal(0);
            double deltaChange = deltaKE - cDynDebugDeltaEnergy(0.0);
            cDynPrintf("ENERGY %12.9f loss= %12.9f TE1= %5.9f TE2= %5.9f PE1= %5.9f KE1= %5.9f PE2= %5.9f KE2= %5.9f\n",
                deltaChange,
                deltaKE,
                cDynDebugPE(0)+cDynDebugKE(0),
                cDynDebugPE(1)+cDynDebugKE(1),
                cDynDebugPE(0), cDynDebugKE(0),
                cDynDebugPE(1), cDynDebugKE(1)
            );
#endif
        }
    }

    // find new world time that is the min of all groups
    cDynTime time=tc;
    for (cDynBaseNode* b=firstBase_;b != NULL; b=b->next_) 
    {
        if (b->status() == CDYN_ACTIVE && b->contact()->time() < time)
            time=b->contact()->time();
    }

#if 0
    if (flag) 
    {
        cDynPrintf("After\n");
        for (cDynBaseNode* b=firstBase_;b != NULL; b=b->next_) 
        {
            b->rootNode()->displaySmallTree(NULL);
        }
    }
#endif

    time_=time;

    // advance f-curves
    if (time-backupLimit_ > backupTime_) 
    {
        backupTime_=time-backupLimit_;
        for (cDynBaseNode* b=firstBase_;b != NULL; b=b->next_) 
        {
            if (b->status() != CDYN_SLEEP && b->status() != CDYN_INACTIVE)
                b->advance(backupTime_);
        }
    }

    // see if any invalid pairs are ok now
    collision_.checkInvalids();

    // update the collision event
    collisionEvent_->time(time+collisionEvent_->inc());

    // process any collision callbacks that are required
    for (cDynBaseNode* b=firstBase_;b != NULL; b=b->next_) 
    {
        b->contact()->event();
    }
    // make call to callback if required
    if (callback_ != NULL) (*callback_)(this);
}

//---------------------------------------------------------------------------
#if 0
//---------------------------------------------------------------------------

size_t cDynWorld::save(void* ptr, size_t size)
{
        char* p=(char *)ptr;
    size_t o= 5*sizeof(cDynTime);
    if (o <= size) {
        *(cDynTime *)p = time_; p += sizeof(cDynTime);
        *(cDynTime *)p = prevTime_; p += sizeof(cDynTime);
        *(cDynTime *)p = backupTime_; p += sizeof(cDynTime);
        *(cDynTime *)p = collisionEvent_->time(); p += sizeof(cDynTime);
        *(cDynTime *)p = collisionEvent_->increment(); p += sizeof(cDynTime);
    }
        size_t s=o;
    p=(char *)ptr;
        for (cDynBaseNode* b=firstBase_;b != NULL; b=b->next_)
                s += b->save((size < s)?p:(p+s), (size < s)?0:(size - s));
        return(s+o);
}

//---------------------------------------------------------------------------

size_t cDynWorld::restore(void* ptr, size_t size)
{
    char* p=(char *)ptr;
    size_t o= 5*sizeof(cDynTime);
    cDynTime w = *(cDynTime *)p; p += sizeof(cDynTime);
    cDynTime prev = *(cDynTime *)p; p += sizeof(cDynTime);
    cDynTime b = *(cDynTime *)p; p += sizeof(cDynTime);
    cDynTime c = *(cDynTime *)p; p += sizeof(cDynTime);
    cDynTime inc = *(cDynTime *)p; p += sizeof(cDynTime);

    cDynTime offset=time_ - w;

    prevTime_ = prev + offset;
    backupTime_ = b + offset;
    collisionEvent_->time(c+offset);
    collisionEvent_->increment(inc);

    size_t s=o;
    p=(char *)ptr;

    cDynBaseNode* base;
    for (base = firstBase_; base != NULL; base = base->next_)
        s += base->restore((size < s) ? p : (p + s), (size < s) ? 0 : (size - s));
    return(s+o);
}

//---------------------------------------------------------------------------

int cDynWorld::check()
{
    int err=0;
    for (cDynBaseNode* b=firstBase_;b != NULL; b=b->next_) 
    {
        err += b->check();
    }
    return(err);
}

//---------------------------------------------------------------------------

bool cDynWorld::checkState()
{
    for (cDynBaseNode* b=firstBase_;b != NULL; b=b->next_) 
    {
        if (!b->checkState()) return(false);
    }
    return(true);
}

//---------------------------------------------------------------------------
#endif
//---------------------------------------------------------------------------
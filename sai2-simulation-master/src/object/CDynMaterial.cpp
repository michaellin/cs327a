//===========================================================================
/*
    This file is part of the dynamics library.
    Copyright (C) 2014, Artificial Intelligence Laboratory,
    Stanford University. All rights reserved.

    \version   $MAJOR.$MINOR.$RELEASE $Rev: 1242 $
*/
//===========================================================================

//---------------------------------------------------------------------------
#include "object/CDynMaterial.h"
#include "global/CDynGlobalDefn.h"
#include "utility/CDynLogger.h"
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
// Known issues:
// the deletion of a material does not necessarily remove all binds from memory.
// Should not cause a problem, although hashes of other materials will be larger
// then strictly required.
//---------------------------------------------------------------------------

cDynMaterialModel cDynMaterial::model_=CDYN_MODEL_MEAN;

int cDynMaterial::h(const cDynMaterial* m)
{
    return(((unsigned long)m/sizeof(cDynMaterial))%size_);
}

//---------------------------------------------------------------------------

cDynMaterial* cDynMaterial::lookup(const cDynMaterial* m)
{
    if (hash_ == NULL) return(NULL);
    
    unsigned int hash=h(m);
    cDynMaterialPair* p = hash_[hash]; 
    while (p != NULL) 
    {
        if (p->a_ == m) return(p->r_);
        p=p->next_;
    }
    
    return(NULL);
}

//---------------------------------------------------------------------------

cDynMaterial* cDynMaterial::create(cDynMaterial* m)
{
    if (hash_ == NULL || m->hash_ == NULL) return(NULL);
    
    unsigned int hash=h(m);
    cDynMaterial* n=new cDynMaterial(0);
    cDynMaterialPair* p = new cDynMaterialPair;
    p->a_=m;
    p->r_=n;
    p->next_=hash_[hash];
    hash_[hash]=p; 

    return(n);
}

//---------------------------------------------------------------------------

void cDynMaterial::unbind(cDynMaterial* m)
{
    cDynMaterial* a;
    cDynMaterial* b;

    // hash to one of the material structures (need better hash)
    if ((unsigned long)this<(unsigned long)m) 
    {
        a=this;b=m;
    } 
    else 
    {
        b=this;a=m;
    }
    
    if (a->hash_ == NULL) return;
    
    unsigned int hash=a->h(b);
    cDynMaterialPair* o=NULL;
    cDynMaterialPair* p = a->hash_[hash]; 
    while (p != NULL) 
    {
        if (p->a_ == b) 
        {
            if (o == NULL) a->hash_[hash]=p->next_;
            else o->next_=p->next_;
            delete p->r_;
            delete p;
            return;
        }
        o=p;
        p=p->next_;
    }
}

//---------------------------------------------------------------------------

cDynMaterial::cDynMaterial(const int size, char* data)
{
    data_=data;
    epsilon_=0.0;
    vlimit_=(double)1e-6;
    for (int i=0;i<CDYN_MATERIAL_FRICTION_MAX;i++)
        friction_[i]=0.0;

    size_=size;
    if (size > 0) 
    {
        hash_= new cDynMaterialPair*[size];
//	    if (cDynDebugLevel == -100)
//		  cDynPrintf("cDynMaterial\n");
//      hash_= &(new (cDynMaterialPair*)[size]);
        for (int i=0;i<size;i++) hash_[i]=NULL;
    } 
    else 
    {
        hash_= NULL;
    }
}

//---------------------------------------------------------------------------

cDynMaterial::~cDynMaterial()
{
    if (size_ > 0) 
    {
        for (int i=0;i< size_;i++) 
        {
            cDynMaterialPair* p=hash_[i];
            while (p != NULL) 
            {
                cDynMaterialPair* t=p;
                p=p->next_;
                delete t->r_;
                delete t;
            }
        }
        delete [] hash_;
    }
}

//---------------------------------------------------------------------------

void cDynMaterial::epsilon(const double e,cDynMaterial* m)
{
    if (m == NULL) { epsilon_=e; return; }
    cDynMaterial* a;
    cDynMaterial* b;

    // hash to one of the material structures (need better hash)
    if ((unsigned long)this<(unsigned long)m) 
    {
        a=this;b=m;
    } 
    else 
    {
        b=this;a=m;
    }

    cDynMaterial* p=a->lookup(b);
    if (p == NULL) p=a->create(b);
    if (p != NULL) p->epsilon_=e;
}

//---------------------------------------------------------------------------

double cDynMaterial::epsilon(cDynMaterial* m)
{
    if (m == NULL) return(epsilon_);
    cDynMaterial* a;
    cDynMaterial* b;

    // hash to one of the material structures (need better hash)
    if ((unsigned long)this<(unsigned long)m) 
    {
        a=this;b=m;
    } 
    else 
    {
        b=this;a=m;
    }
    
    cDynMaterial* p=a->lookup(b);
    if (p != NULL) 
    {
        return(p->epsilon());
    } 
    else 
    {
        double ea=a->epsilon();
        double eb=b->epsilon();
        switch (model_) 
        {
            case CDYN_MODEL_MIN: return((ea<eb)?ea:eb);
            case CDYN_MODEL_MEAN: return((ea+eb)*0.5f);
            case CDYN_MODEL_MAX: return((ea<eb)?eb:ea);
        }
    }

    return(0.0); 
}

//---------------------------------------------------------------------------

void cDynMaterial::friction(const cDynFrictionType type, const double u,cDynMaterial* m)
{
    if (m == NULL) { friction_[type]=u; return; }

    cDynMaterial* a;
    cDynMaterial* b;

    // hash to one of the material structures (need better hash)
    if ((unsigned long)this<(unsigned long)m) 
    {
        a=this;b=m;
    } 
    else 
    {
        b=this;a=m;
    }
    
    cDynMaterial* p=a->lookup(b);
    if (p == NULL) p=a->create(b);
    if (p != NULL) p->friction_[type]=u;
}

//---------------------------------------------------------------------------

double cDynMaterial::friction(const cDynFrictionType type, cDynMaterial* m)
{
    if (m == NULL) return(friction_[type]);
    cDynMaterial* a;
    cDynMaterial* b;

    // hash to one of the material structures (need better hash)
    if ((unsigned long)this<(unsigned long)m) 
    {
        a=this;b=m;
    } 
    else 
    {
        b=this;a=m;
    }

    cDynMaterial* p=a->lookup(b);
    if (p != NULL) 
    {
        return(p->friction(type));
    } 
    else 
    {
        double ua=a->friction(type);
        double ub=b->friction(type);
        switch (model_) 
        {
            case CDYN_MODEL_MIN: return((ua<ub)?ua:ub);
            case CDYN_MODEL_MEAN: return((ua+ub)*0.5f);
            case CDYN_MODEL_MAX: return((ua<ub)?ub:ua);
        }
    }

    return(0.0f); 
}

//---------------------------------------------------------------------------

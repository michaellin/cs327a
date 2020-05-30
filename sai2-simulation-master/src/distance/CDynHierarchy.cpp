//===========================================================================
/*
    This file is part of the dynamics library.
    Copyright (C) 2014, Artificial Intelligence Laboratory,
    Stanford University. All rights reserved.

    \version   $MAJOR.$MINOR.$RELEASE $Rev: 1242 $
*/
//===========================================================================

//---------------------------------------------------------------------------
#include "distance/CDynHierarchy.h"
#include "utility/CDynLogger.h"
#include "distance/CDynBSphere.h"
#include <stdlib.h>
//---------------------------------------------------------------------------
static void _cDynHierarchyCheck(cDynBSphere* a, cDynBSphere* b);
//---------------------------------------------------------------------------

static void _cDynHierarchyCheck(cDynBSphere* a, cDynBSphere* b)
{
    if (CDYN_BSPHERE_LEFT(b) == NULL && CDYN_BSPHERE_RIGHT(b) == NULL) 
    {
        cDynVector3 d;
        d.subtract(CDYN_BSPHERE_P(b), CDYN_BSPHERE_P(a));
        double dc=d.magnitude() + CDYN_BSPHERE_R(b);
        if (dc > CDYN_BSPHERE_R(a) + 1e-8) 
        {
            cDynPrintf("Error Ancestor does not contain all children\n");
            exit(1);
        }
    } 
    else 
    {
        if (CDYN_BSPHERE_LEFT(b) != NULL) _cDynHierarchyCheck(a,CDYN_BSPHERE_LEFT(b));
        if (CDYN_BSPHERE_RIGHT(b) != NULL) _cDynHierarchyCheck(a,CDYN_BSPHERE_RIGHT(b));
    }
}

//---------------------------------------------------------------------------

cDynBSphere *cDynHierarchy(cDynBSphere *list)
{
    static int call=0;

    call++;
    if (list == NULL || CDYN_BSPHERE_LEFT(list) == NULL) return(list);

    cDynVector3 avg=CDYN_BSPHERE_P(list);
    cDynVector3 min;
    min.add(CDYN_BSPHERE_P(list), -CDYN_BSPHERE_R(list));
    cDynVector3 max;
    max.add(CDYN_BSPHERE_P(list), CDYN_BSPHERE_R(list));
    int n=1;

    cDynVector3 tmpV;
    for (cDynBSphere* s=CDYN_BSPHERE_LEFT(list); s != NULL; s=CDYN_BSPHERE_LEFT(s)) 
    {
        avg += CDYN_BSPHERE_P(s);
        tmpV.add(CDYN_BSPHERE_P(s), -CDYN_BSPHERE_R(s));
        min.minimum(tmpV);
        tmpV.add(CDYN_BSPHERE_P(s), -CDYN_BSPHERE_R(s));
        max.maximum(tmpV);
        n++;
    }

    avg *= (1.0f/(double)n);
    cDynVector3 center;
    center.add(min, max);
    center *= 0.5f;

    // heuristic 2	
    double r2=0.0f;
    for (cDynBSphere* s=list; s != NULL; s=CDYN_BSPHERE_LEFT(s)) 
    {
        cDynVector3 d;
        d.subtract(CDYN_BSPHERE_P(s), center);
        double r = d.magnitude() + CDYN_BSPHERE_R(s);
        if (r > r2) r2=r;
    }

    cDynVector3 range;
    range.subtract(max, min);

    // find longest axis
    int c;
    if (range[0] > range[1] && range[0] > range[2]) 
        c=0; // x range largest
    else if (range[1] > range[2])			
        c=1; // y range largest
    else						
        c=2; // z range largest

    // split spheres along mean of longest axes
    int toggle=0;
    cDynBSphere *listl=NULL;
    cDynBSphere *listr=NULL;
    
    for (cDynBSphere* s=list;s != NULL;) 
    {
        cDynBSphere* snext=CDYN_BSPHERE_LEFT(s);
        if (CDYN_BSPHERE_P(s)[c] == avg[c]) 
        {
            if (toggle>0) 
            {
                CDYN_BSPHERE_LEFT(s)=listl; listl=s; toggle--;
            } 
            else 
            {
                CDYN_BSPHERE_LEFT(s)=listr; listr=s; toggle++;
            }
        } 
        else if (CDYN_BSPHERE_P(s)[c] < avg[c]) 
        {
            CDYN_BSPHERE_LEFT(s)=listl;
            listl=s;
            toggle--;
        } 
        else 
        {
            CDYN_BSPHERE_LEFT(s)=listr;
            listr=s;
            toggle++;
        }
        s=snext;
    }

    // fixup if odd round off errors occur
    if (listl == NULL) 
    {
        listl=listr; listr=CDYN_BSPHERE_LEFT(listr); CDYN_BSPHERE_LEFT(listl)=NULL;
    } 
    else if (listr == NULL) 
    {
        listr=listl; listl=CDYN_BSPHERE_LEFT(listl); CDYN_BSPHERE_LEFT(listr)=NULL;
    }

    // recurse left and right lists
    cDynBSphere* left=cDynHierarchy(listl);
    cDynBSphere* right=cDynHierarchy(listr);

    // heuristic 1
    cDynVector3 d;
    d.subtract(CDYN_BSPHERE_P(right), CDYN_BSPHERE_P(left));
    double dis=d.magnitude();
    double r1 = 0.5f*(dis
        + (( CDYN_BSPHERE_R(left) > CDYN_BSPHERE_R(right) - dis)?CDYN_BSPHERE_R(left):CDYN_BSPHERE_R(right) - dis)
        + ((CDYN_BSPHERE_R(right) > CDYN_BSPHERE_R(left) - dis )?CDYN_BSPHERE_R(right):CDYN_BSPHERE_R(left) - dis));

    // use best heuristic to compute p and r for new sphere
    cDynVector3 p;
    double r;
    if (r1 < r2) 
    {
        if (dis < 1e-4) 
        { 
            // if spheres very nearly overlap
            r=r1+0.5f*dis;
            p=CDYN_BSPHERE_P(left);
        } 
        else 
        {
            d *= (1.0f/dis);
            r = r1;
            double s = r - (( CDYN_BSPHERE_R(left) > CDYN_BSPHERE_R(right) - dis)?CDYN_BSPHERE_R(left) :CDYN_BSPHERE_R(right) - dis);
            p.multiply(d, s);
            p+=CDYN_BSPHERE_P(left);
        }
    } 
    else 
    {
        p=center;
        r=r2;
    }

    cDynBSphere* node=new cDynBSphere(&p,r,NULL);
    CDYN_BSPHERE_LEFT(node)=left;
    CDYN_BSPHERE_RIGHT(node)=right;

#if 0
    // promote prim pointer if both left and right branches all point to same prim
    if (CDYN_BSPHERE_LEFT(node) != NULL && CDYN_BSPHERE_RIGHT(node) != NULL
     && node->CDYN_BSPHERE_PRIM(left) != NULL && node->CDYN_BSPHERE_PRIM(left) == node->CDYN_BSPHERE_PRIM(right)) 
    {
        CDYN_BSPHERE_PRIM(node) = node->CDYN_BSPHERE_PRIM(left);
    }
#endif

    if (r - CDYN_BSPHERE_R(left) < 1e-4) 
    {
        CDYN_BSPHERE_STATE(left)=0;
    }
    if (r - CDYN_BSPHERE_R(right) < 1e-4) 
    {
        CDYN_BSPHERE_STATE(right)=0;
    }

    return(node);
}

//---------------------------------------------------------------------------

void cDynHierarchyFullCheck(cDynBSphere* root)
{
    if (CDYN_BSPHERE_LEFT(root) != NULL) _cDynHierarchyCheck(root,CDYN_BSPHERE_LEFT(root));
    if (CDYN_BSPHERE_RIGHT(root) != NULL) _cDynHierarchyCheck(root,CDYN_BSPHERE_RIGHT(root));
    if (CDYN_BSPHERE_LEFT(root) != NULL) cDynHierarchyFullCheck(CDYN_BSPHERE_LEFT(root));
    if (CDYN_BSPHERE_RIGHT(root) != NULL) cDynHierarchyFullCheck(CDYN_BSPHERE_RIGHT(root));
}

//---------------------------------------------------------------------------

int cDynHierarchyLeafCount(cDynBSphere* root)
{
    if (root == NULL) return(0);
    if (CDYN_BSPHERE_LEFT(root) == NULL && CDYN_BSPHERE_RIGHT(root) == NULL) return(1);
    return(cDynHierarchyLeafCount(CDYN_BSPHERE_LEFT(root)) + cDynHierarchyLeafCount(CDYN_BSPHERE_RIGHT(root)));
}

//---------------------------------------------------------------------------

int cDynHierarchyCount(cDynBSphere* root)
{
    if (root == NULL) return(0);
    if (CDYN_BSPHERE_LEFT(root) == NULL && CDYN_BSPHERE_RIGHT(root) == NULL) return(1);
    return(1 + cDynHierarchyCount(CDYN_BSPHERE_LEFT(root)) + cDynHierarchyCount(CDYN_BSPHERE_RIGHT(root)));
}

//---------------------------------------------------------------------------
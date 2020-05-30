//===========================================================================
/*
    This file is part of the dynamics library.
    Copyright (C) 2014, Artificial Intelligence Laboratory,
    Stanford University. All rights reserved.

    \version   $MAJOR.$MINOR.$RELEASE $Rev: 1242 $
*/
//===========================================================================

//---------------------------------------------------------------------------
#include <stdio.h>
//---------------------------------------------------------------------------
#include "extras/CDynHull.h"
#include "object/CDynObject.h"
#include "distance/CDynPrim.h"
#include "object/CDynObject.h"
//---------------------------------------------------------------------------
extern "C" 
{
#include "qhull_a.h"
}
//---------------------------------------------------------------------------

static int _cDynPlaner(coordT **rows,int n, coordT err)
{
    coordT mean[3];
    coordT normal[3];
    {
        for (int i=0;i<3;i++) 
        {
            normal[i]=mean[i]=0.0;
        }
    }
    
    {
        for (int i = 0; i < n; i++) 
        {
            int j;
            coordT nv[3];

            if (i == n-1) 
                j = 0;
            else 
                j = i+1;

            nv[0] = (rows[i][1] - rows[j][1]) * (rows[i][2] + rows[j][2]);
            nv[1] = (rows[i][2] - rows[j][2]) * (rows[i][0] + rows[j][0]);
            nv[2] = (rows[i][0] - rows[j][0]) * (rows[i][1] + rows[j][1]);

            coordT dot=0.0;
            for (int k=0;k<3;k++) 
            {
                dot += nv[k]*normal[k];
                mean[k] += rows[i][k];
            }
            if (dot < 0.0)
            {
                for (int k=0;k<3;k++) 
                    normal[k] -= nv[k];
            }
            else
            {
                for (int k=0;k<3;k++) 
                    normal[k] += nv[k];
            }
        }
    }

    {
        for (int k=0;k<3;k++)
            mean[k] /= n;
    }

    coordT d=0.0;
    {
        for (int k=0;k<3;k++) 
            d += mean[k]*normal[k];
    }

    {
        for (int i = 0; i < n; i++) 
        {
            if (fabs(rows[i][0]*normal[0] + rows[i][1]*normal[1] + rows[i][2]*normal[2] - d) > err) 
                return(-1);
        }
    }

    int max=0;
    {
        for (int k=1;k<3;k++)
            if (fabs(normal[k]) > fabs(normal[max])) 
                max=k;
    }

    return(max);
}

//---------------------------------------------------------------------------

void cDynHull(cDynObject* obj)
{
    int n=0;
    int np=0;
    double r=0.0;
    double err=0.0;
    cDynMaterial* mat;
    const cDynPrim* p;

    for (p=obj->geometry.prim();p != NULL; p=p->next) 
    {
        if (p->r > r) 
        {
            r=p->r;
            err=p->err;
            mat=p->m;
        }
        n += p->num;
        np++;
    }

    if (n <= 3) 
    {
        cDynPrimitive* h=obj->geometry.begin(CDYN_HULL);
        h->radius(r);
        h->error(err);
        h->material(mat);
        for (const cDynPrim* p=obj->geometry.prim();p != NULL; p=p->next) 
        {
            for (int i=0;i<p->num;i++)
                h->vertex(p->v[i]);
        }
        h->end();
        obj->geometry.end();
        return;
    }
    coordT *points= new coordT[3*n];
    coordT **rows= new coordT*[n];
    vertexT *vertex;
    int curlong, totlong;

    int k=0;
    for (p=obj->geometry.prim();p != NULL; p=p->next) 
    {
        for (int i=0;i<p->num;i++)
            for (int j=0;j<3;j++)
                points[k++]=p->v[i][j];
    }
    
    for (int i=0; i < n; i++)
    {	
        rows[i]= points+3*i;
    }

    char options[255];
    int plane= _cDynPlaner(rows,n,err*0.25);
    if (plane == -1) 
    {
        sprintf(options,"qhull p QJ C-%0.15f Pp",err);
    } 
    else 
    {
        sprintf(options,"qhull p QJ C-%0.15f Qb%d:0B%d:0 Pp",err,plane,plane);
    }
    
    if (!qh_new_qhull (3, n, points, False, options, NULL, NULL)) 
    {
        cDynGeometry geom;
        cDynPrimitive* h= geom.begin(CDYN_HULL);
        h->radius(r);
        h->error(err);
        h->material(mat);
        FORALLvertices 
        {
            h->vertex(
                rows[qh_pointid(vertex->point)][0],
                rows[qh_pointid(vertex->point)][1],
                rows[qh_pointid(vertex->point)][2]
            );
        }
        h->end();
        geom.end();
        obj->geometry = geom;
    }

    qh_freeqhull(!qh_ALL);                    
    qh_memfreeshort (&curlong, &totlong);    
    
    delete [] points;
    delete [] rows;
}

//---------------------------------------------------------------------------


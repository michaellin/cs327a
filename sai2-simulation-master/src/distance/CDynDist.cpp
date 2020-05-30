//===========================================================================
/*
    This file is part of the dynamics library.
    Copyright (C) 2014, Artificial Intelligence Laboratory,
    Stanford University. All rights reserved.

    \version   $MAJOR.$MINOR.$RELEASE $Rev: 1242 $
*/
//===========================================================================

//---------------------------------------------------------------------------
#include "distance/CDynDist.h"
#include "matrix/CDynVector3.h"
#include "utility/CDynLogger.h"
//---------------------------------------------------------------------------
#define CDYN_DIST_DOT(v1,v2)			cDotV3V3(v1,v2)
#define CDYN_DIST_MAGNITUDE(v)			cMagnitudeV3(v)
#define CDYN_DIST_COPY(res,v)			cSetV3V3(res,v)
#define CDYN_DIST_COPY_NEG(res,v)		cNegateV3V3(res,v)
#define CDYN_DIST_SUBTRACT(res,v1,v2)   cSubV3V3V3(res,v1,v2)
#define CDYN_DIST_MULTIPLYEQUAL(res,c)	cMulV3S1(res,c)
//---------------------------------------------------------------------------
// internal variables
static double _cDynDistEps = CDYN_DIST_EPS;
static int _cDynDistTnum = 0;
//---------------------------------------------------------------------------
// internal functions
//static int _cDynDistCSFCN(const cDynVector3 *vec,int num,const double eta[3],double *suf);
//---------------------------------------------------------------------------

static double _cDynDistDSBP(int *nov,int v_index[2][4],double lambda[4],
                               double y[4][3],double del[4][4],
                               double zsol[3],bool backup);

static double _cDynDistDSBP2(int *nov,int v_index[2][4],double lambda[4],
                               double y[4][3],double del[4][4],
                               double zsol[3]);

//---------------------------------------------------------------------------

void cDynDistNormal(const cDynVector3 *veca, 
    const cDynVector3 *vecb, 
    const cDynDistWitness *dw, 
    cDynVector3 *n)
{
    int i;
    double d,lambda;
    double *f = (double *)n;
    f[0]=f[1]=f[2]=0;
    for(i=0;i<CDYN_DISTW_NOV(dw);i++)
    {
        const double *va = veca[CDYN_DISTW_INDEX(dw)[0][i]];
        const double *vb = vecb[CDYN_DISTW_INDEX(dw)[1][i]];
        lambda = CDYN_DISTW_LAMBDA(dw)[i];
        f[0] += (va[0]-vb[0])*lambda;
        f[1] += (va[1]-vb[1])*lambda;
        f[2] += (va[2]-vb[2])*lambda;
    }
    d = CDYN_DIST_MAGNITUDE(f);
    if (d > CDYN_DIST_EPS)
        CDYN_DIST_MULTIPLYEQUAL(f,1/d);
}

//---------------------------------------------------------------------------

void cDynDistPoint(const cDynVector3 *vec, 
    const cDynDistWitness *dw, 
    const int index, 
    cDynVector3 *x)
{
    int i;
    double lambda;
    double *f = (*x);
    f[0]=f[1]=f[2]=0;
    for(i=0;i<CDYN_DISTW_NOV(dw);i++)
    {
        const double *v = vec[CDYN_DISTW_INDEX(dw)[index][i]];
        lambda = CDYN_DISTW_LAMBDA(dw)[i];
        f[0] += v[0]*lambda;
        f[1] += v[1]*lambda;
        f[2] += v[2]*lambda;
    }
}

//---------------------------------------------------------------------------

static int _cDynDistCSFCN(const cDynVector3 *vec, 
    int num, 
    double eta[3], 
    double *suf)
{
    int i;
    int index =0;
    double max;
    *suf = CDYN_DIST_DOT(vec[0],eta);
    for(i=1;i<num;i++)
    {	
        max = CDYN_DIST_DOT(vec[i],eta);
        if (max > *suf)
        {
            *suf = max;
            index = i;
        }
    }
    return index;
}

//---------------------------------------------------------------------------

double cDynDistDistance2(const cDynVector3* vA,
    const int na,
    const cDynVector3* vB,
    const int nb,cDynDistWitness *dw,
    const double max,
    const bool useMax)
{
    int i,j,ri,rj;
    double lastdstsq, dstsq, si, sj, g, st, dst;
    double del[4][4];
    double yf[4][3];
    double zsolf[3],nzsolf[3];

    if (CDYN_DISTW_NOV(dw)==0) 
    {
        CDYN_DISTW_NOV(dw)=1;
        CDYN_DISTW_INDEX(dw)[0][0]=0;
        CDYN_DISTW_INDEX(dw)[1][0]=0;
        CDYN_DISTW_LAMBDA(dw)[0]=1;
    }

    // update polytope and dot product table
    for (i=0;i<CDYN_DISTW_NOV(dw);i++)
        CDYN_DIST_SUBTRACT(yf[i],vA[CDYN_DISTW_INDEX(dw)[0][i]],vB[CDYN_DISTW_INDEX(dw)[1][i]]);

    for (i=0;i<CDYN_DISTW_NOV(dw);i++)
        for (j=0;j<=i;j++)
            del[i][j] = CDYN_DIST_DOT(yf[i],yf[j]);

    lastdstsq = del[0][0]+del[0][0]+1;

    while (true) 
    {
        // call distance sub algorithm : backup is always false
        dstsq = _cDynDistDSBP2(&(CDYN_DISTW_NOV(dw)),CDYN_DISTW_INDEX(dw),CDYN_DISTW_LAMBDA(dw),yf,del,zsolf);

        if (CDYN_DISTW_NOV(dw) == 4)
        {
            // penetration --> no backup
            return 0;
        }

        if (dstsq >= lastdstsq)
        {
            // no more progress or some kind of error: just return --> no backup
            return cDynSqrt(dstsq);
        }

        lastdstsq=dstsq;

        CDYN_DIST_COPY_NEG(nzsolf,zsolf);  // nzsolf = - zsolf

        ri = _cDynDistCSFCN(vA,na,nzsolf,&si);
        rj = _cDynDistCSFCN(vB,nb,zsolf,&sj);

        g = dstsq + si + sj;

        // see if we are close enough to the solution
        if (g <= _cDynDistEps) 
        {
            // found solution, return
            return cDynSqrt(dstsq);
        }

        if (useMax)
        {
            dst = cDynSqrt(dstsq);
            st = -(si+sj)/dst;
            if (st > max) 
            {
                // found penetration into the max boundary???
                return(st);
            }
        }

        // if not add new point to polytope and try again
        // move first point to last position
        CDYN_DISTW_INDEX(dw)[0][CDYN_DISTW_NOV(dw)]=CDYN_DISTW_INDEX(dw)[0][0];
        CDYN_DISTW_INDEX(dw)[1][CDYN_DISTW_NOV(dw)]=CDYN_DISTW_INDEX(dw)[1][0];
        CDYN_DIST_COPY(yf[CDYN_DISTW_NOV(dw)],yf[0]);
        for (i=1;i<CDYN_DISTW_NOV(dw);i++) 
            del[CDYN_DISTW_NOV(dw)][i] = del[i][0];
        del[CDYN_DISTW_NOV(dw)][CDYN_DISTW_NOV(dw)]=del[0][0];

        // put new point in first spot
        CDYN_DISTW_INDEX(dw)[0][0]=ri;
        CDYN_DISTW_INDEX(dw)[1][0]=rj;

        CDYN_DIST_SUBTRACT(yf[0],vA[ri],vB[rj]);

        // update dot product table
        for (i=0;i<=CDYN_DISTW_NOV(dw);i++)
            del[i][0] = CDYN_DIST_DOT(yf[i],yf[0]);
        CDYN_DISTW_NOV(dw)++;
    }
    // penetration --> this line can't be reached.
    return 0;
}

//---------------------------------------------------------------------------

double cDynDistDistance(const cDynVector3* vA,
    const int na,
    const cDynVector3* vB,
    const int nb,cDynDistWitness *dw,
    const double max,
    const bool useMax)
{
    bool backup=false;
    int iter=0;
    int oldnov=0;
    int ncy=0;
    int ri,rj,i,j,k,l,ii,kk,ll;
    int oldv_index[2][4], iord[4];
    double lastdstsq, dstsq, si, sj, g;
    double del[4][4],olddel[4][4];
    double yf[4][3],oldyf[4][3];
    double zsolf[3],nzsolf[3];

    _cDynDistTnum++;

    if (CDYN_DISTW_NOV(dw)==0) 
    {
        CDYN_DISTW_NOV(dw)=1;
        CDYN_DISTW_INDEX(dw)[0][0]=0;
        CDYN_DISTW_INDEX(dw)[1][0]=0;
        CDYN_DISTW_LAMBDA(dw)[0]=1;
    }

    // update polytope and dot product table
    if (useMax)
    {
        for (i=0;i<CDYN_DISTW_NOV(dw);i++)
            //objA->v[CDYN_DISTW_INDEX(dw)[0][i]].subtract(objB->v[CDYN_DISTW_INDEX(dw)[1][i]],y[i]);
            CDYN_DIST_SUBTRACT(yf[i],vA[CDYN_DISTW_INDEX(dw)[0][i]],vB[CDYN_DISTW_INDEX(dw)[1][i]]);
    }
    else
    {
        for (i=0;i<CDYN_DISTW_NOV(dw);i++)
            //objA->v[i].subtract(objB->v[i],y[i]);
            CDYN_DIST_SUBTRACT(yf[i],vA[i],vB[i]);
    }

    for (i=0;i<CDYN_DISTW_NOV(dw);i++)
        for (j=0;j<=i;j++)
            del[i][j] = CDYN_DIST_DOT(yf[i],yf[j]);

    lastdstsq = del[0][0]+del[0][0]+1;

    while (true) 
    {
        iter++;
        ncy++;

        // call distance sub algorithm
        dstsq = _cDynDistDSBP(&(CDYN_DISTW_NOV(dw)),CDYN_DISTW_INDEX(dw),CDYN_DISTW_LAMBDA(dw),yf,del,zsolf,backup);

        if (dstsq >= lastdstsq || CDYN_DISTW_NOV(dw) == 4) 
        {
            if (backup) 
            {
                return cDynSqrt(dstsq);
            }
            backup=true;
            if (ncy == 1) 
                continue;
            //		    exit(-1);
            CDYN_DISTW_NOV(dw) = oldnov;
            for (k=0;k<CDYN_DISTW_NOV(dw);k++) 
            {
                CDYN_DISTW_INDEX(dw)[0][k] = oldv_index[0][k];
                CDYN_DISTW_INDEX(dw)[1][k] = oldv_index[1][k];
                CDYN_DIST_COPY(yf[k],oldyf[k]);

                for (l=0;l<=k;l++)
                    del[k][l]=olddel[k][l];
            }
            //memcpy(del, olddel, sizeof(del));
            continue;
        }
        lastdstsq=dstsq;

        // find new point to add to polytope
        // nzsolf = -zsolf
        CDYN_DIST_COPY_NEG(nzsolf,zsolf);

        ri = _cDynDistCSFCN(vA,na,nzsolf,&si);
        rj = _cDynDistCSFCN(vB,nb,zsolf,&sj);

        g = dstsq + si + sj;

        // see if we are close enough to the solution
        if (g <= _cDynDistEps) 
        {
            // found solution, return
            return cDynSqrt(dstsq);
        }

        if (useMax)
        {
            double dst = cDynSqrt(dstsq);
            double st = -(si+sj)/dst;
            if (st > max) 
                return(st);
        }

        // if not add new point to polytope and try again
        // move first point to last position
        CDYN_DISTW_INDEX(dw)[0][CDYN_DISTW_NOV(dw)]=CDYN_DISTW_INDEX(dw)[0][0];
        CDYN_DISTW_INDEX(dw)[1][CDYN_DISTW_NOV(dw)]=CDYN_DISTW_INDEX(dw)[1][0];
        CDYN_DIST_COPY(yf[CDYN_DISTW_NOV(dw)],yf[0]);
        for (i=1;i<CDYN_DISTW_NOV(dw);i++) 
            del[CDYN_DISTW_NOV(dw)][i] = del[i][0];
        del[CDYN_DISTW_NOV(dw)][CDYN_DISTW_NOV(dw)]=del[0][0];

        // put new point in first spot
        CDYN_DISTW_INDEX(dw)[0][0]=ri;
        CDYN_DISTW_INDEX(dw)[1][0]=rj;

        CDYN_DIST_SUBTRACT(yf[0],vA[ri],vB[rj]);

        // update dot product table
        for (i=0;i<=CDYN_DISTW_NOV(dw);i++)
            del[i][0] = CDYN_DIST_DOT(yf[i],yf[0]);
        CDYN_DISTW_NOV(dw)++;

        // save old values of nov, v_index[0], v_index[1], y and del
        oldnov=CDYN_DISTW_NOV(dw);
        for (k=0;k<CDYN_DISTW_NOV(dw);k++) 
        {
            oldv_index[0][k]=CDYN_DISTW_INDEX(dw)[0][k];
            oldv_index[1][k]=CDYN_DISTW_INDEX(dw)[1][k];
            CDYN_DIST_COPY(oldyf[k],yf[k]);

            for (l=0;l<=k;l++)
                olddel[k][l]=del[k][l];
        }
        //memcpy(olddel, del, CDYN_DISTW_NOV(dw) * 4 * sizeof(double));

        // if CDYN_DISTW_NOV(dw) == 4, 
        // rearrange del[1][0], del[2][0] and del[3][0]
        // in non decreasing order
        if (CDYN_DISTW_NOV(dw) == 4) 
        {
            iord[0] = 0;
            iord[1] = 1;
            iord[2] = 2;
            if (del[2][0] < del[1][0]) 
            {
                iord[1] = 2;
                iord[2] = 1;
            }
            ii = iord[1];
            if (del[3][0] < del[ii][0]) 
            {
                iord[3] = iord[2];
                iord[2] = iord[1];
                iord[1] = 3;
            } 
            else 
            {
                ii = iord[2];
                if (del[3][0] < del[ii][0]) 
                {
                    iord[3] = iord[2];
                    iord[2] = 3;
                } 
                else 
                {
                    iord[3] = 3;
                }
            }
            // reorder v_index[0],v_index[1] y and del
            for (k=1;k<CDYN_DISTW_NOV(dw);k++) 
            {
                kk = iord[k];
                CDYN_DISTW_INDEX(dw)[0][k] = oldv_index[0][kk];
                CDYN_DISTW_INDEX(dw)[1][k] = oldv_index[1][kk];
                CDYN_DIST_COPY(yf[k],oldyf[kk]);
                for (l=0;l<k;l++) 
                {
                    ll = iord[l];
                    if (kk >= ll)
                        del[k][l] = olddel[kk][ll];
                    else
                        del[k][l] = olddel[ll][kk];
                }
                del[k][k] = olddel[kk][kk];
            }
        }
    }
  return 0;
}

//---------------------------------------------------------------------------

static double _cDynDistDSBP(int *nov,int v_index[2][4],double lambda[4],
                           double y[4][3],double del[4][4],
                           double zsol[3],bool backup)
{
    int k,l,kk,ll;
    int novd=0;
    int v_indexd[2][4],iord[4];
    double sum;
    double e132,e142,e123,e143,e213,e243;
    double e124,e134,e214,e234,e314,e324;
    double d1[15],d2[15],d3[15],d4[15];
    double yd[4][3]; 
    double deld[4][4];
    double zsold[3];
    double dstsq=0.0f;
    double alsd[4], dstsqd;


    d1[0]=d2[1]=d3[3]=d4[7]=1.0f;
    if (!backup) 
    {
    // regular distance subalgoritm begins.
    switch (*nov) 
    {
        case 1: 
        {       
            //
            // case of a single point.
            //
            lambda[0] = d1[0];
            CDYN_DIST_COPY(zsol,y[0]);
            dstsq = del[0][0];
            return(dstsq);
            break;
            // END case single points
        } 
    
        case 2: 
        {
            //
            // case of two points
            //

            // check optimality of vertex 1
            d2[2] = del[0][0] - del[1][0];
            if (d2[2] <= 0.0f) 
            {
                *nov = 1;
                lambda[0] = d1[0];
                CDYN_DIST_COPY(zsol,y[0]);
                dstsq = del[0][0];
                return(dstsq);
            }

            // check optimality of line segment 1-2
            d1[2] = del[1][1] - del[1][0];
            if (d1[2]>0.0f && d2[2]>0.0f) 
            {
                sum = d1[2] + d2[2];
                lambda[0] = d1[2]/sum;
                lambda[1] = 1.0f - lambda[0];
                zsol[0] = y[1][0] + lambda[0]*(y[0][0] - y[1][0]);
                zsol[1] = y[1][1] + lambda[0]*(y[0][1] - y[1][1]);
                zsol[2] = y[1][2] + lambda[0]*(y[0][2] - y[1][2]);
                dstsq = CDYN_DIST_DOT(zsol,zsol);
                return(dstsq);
            }
            
            // check optimality of vertex 2
            if (d1[2]<=0.0f) 
            {
                *nov = 1;
                v_index[0][0] = v_index[0][1];
                v_index[1][0] = v_index[1][1];
                lambda[0] = d2[1];
                CDYN_DIST_COPY(zsol,y[1]);
                dstsq = del[1][1];
                CDYN_DIST_COPY(y[0],y[1]);
                del[0][0] = del[1][1];
                return(dstsq);
            }
            break;
        }  // END case two points

        case 3: 
        {   
            //
            // case of three points
            //

            // check optimality of vertex 1.
            d2[2] = del[0][0] - del[1][0];
            d3[4] = del[0][0] - del[2][0];
            if (d2[2]<=0.0f && d3[4]<=0.0f)
            {
                *nov = 1;
                lambda[0] = d1[0];
                CDYN_DIST_COPY(zsol,y[0]);
                dstsq = del[0][0];
                return(dstsq);
            }
                
            // check optimality of line segment 1-2
            e132 = del[1][0] - del[2][1];
            d1[2] = del[1][1] - del[1][0];
            d3[6] = d1[2]*d3[4] + d2[2]*e132;
            if (d1[2]>0.0f && d2[2]> 0.0f && d3[6]<=0.0f) 
            {
                *nov = 2;
                sum = d1[2] + d2[2];
                lambda[0] = d1[2]/sum;
                lambda[1] = 1.0f - lambda[0];
                zsol[0] = y[1][0] + lambda[0]*(y[0][0] - y[1][0]);
                zsol[1] = y[1][1] + lambda[0]*(y[0][1] - y[1][1]);
                zsol[2] = y[1][2] + lambda[0]*(y[0][2] - y[1][2]);
                dstsq = CDYN_DIST_DOT(zsol,zsol);
                return(dstsq);
            }
            
            // check optimality of line segment 1-3
            e123 = del[2][0] - del[2][1];
            d1[4] = del[2][2] - del[2][0];
            d2[6] = d1[4]*d2[2] + d3[4]*e123;
            if (d1[4]>0.0f && d2[6]<=0.0f && d3[4]>0.0f) 
            {
                *nov = 2;
                v_index[0][1] = v_index[0][2];
                v_index[1][1] = v_index[1][2];
                sum = d1[4] + d3[4];
                lambda[0] = d1[4]/sum;
                lambda[1] = 1.0f - lambda[0];
                zsol[0] = y[2][0] + lambda[0]*(y[0][0] - y[2][0]);
                zsol[1] = y[2][1] + lambda[0]*(y[0][1] - y[2][1]);
                zsol[2] = y[2][2] + lambda[0]*(y[0][2] - y[2][2]);
                dstsq = CDYN_DIST_DOT(zsol,zsol);
                CDYN_DIST_COPY(y[1],y[2]);
                del[1][0] = del[2][0];
                del[1][1] = del[2][2];
                return(dstsq);
            }
                
            // check optimality of face 1-2-3
            e213 = -e123;
            d2[5] = del[2][2] - del[2][1];
            d3[5] = del[1][1] - del[2][1];
            d1[6] = d2[5]*d1[2] + d3[5]*e213;
            if (d1[6]>0.0f && d2[6]>0.0f && d3[6]>0.0f) 
            {
                sum = d1[6] + d2[6] + d3[6];
                lambda[0] = d1[6]/sum;
                lambda[1] = d2[6]/sum;
                lambda[2] = 1.0f - lambda[0] - lambda[1];
                zsol[0] = y[2][0] + lambda[0]*(y[0][0] - y[2][0]) + lambda[1]*(y[1][0] - y[2][0]);
                zsol[1] = y[2][1] + lambda[0]*(y[0][1] - y[2][1]) + lambda[1]*(y[1][1] - y[2][1]);
                zsol[2] = y[2][2] + lambda[0]*(y[0][2] - y[2][2]) + lambda[1]*(y[1][2] - y[2][2]);
                dstsq = CDYN_DIST_DOT(zsol,zsol);
                //	  cDynPrintf("case 3: sum=%13.9f\n",sum);
                return(dstsq);
            }
                
            // check optimality of vertex 2.
            if (d1[2]<=0.0f && d3[5]<=0.0f) 
            {
                *nov = 1;
                v_index[0][0] = v_index[0][1];
                v_index[1][0] = v_index[1][1];
                lambda[0] = d2[1];
                CDYN_DIST_COPY(zsol,y[1]);
                dstsq = del[1][1];
                CDYN_DIST_COPY(y[0],y[1]);
                del[0][0] = del[1][1];
                return(dstsq);
            }

            // check optimality of vertex 3
            if (d1[4]<=0.0f && d2[5]<=0.0f) 
            {
                *nov = 1;
                v_index[0][0] = v_index[0][2];
                v_index[1][0] = v_index[1][2];
                lambda[0] = d3[3];
                CDYN_DIST_COPY(zsol,y[2]);
                dstsq = del[2][2];
                CDYN_DIST_COPY(y[0],y[2]);
                del[0][0] = del[2][2];
                return(dstsq);
            }

            // check optimality of line segment 2-3
            if (d1[6]<=0.0f && d2[5]>0.0f && d3[5]>0.0f) 
            {
                *nov = 2;
                v_index[0][0] = v_index[0][2];
                v_index[1][0] = v_index[1][2];
                sum = d2[5] + d3[5];
                lambda[1] = d2[5]/sum;
                lambda[0] = 1.0f - lambda[1];
                zsol[0] = y[2][0] + lambda[1]*(y[1][0] - y[2][0]);
                zsol[1] = y[2][1] + lambda[1]*(y[1][1] - y[2][1]);
                zsol[2] = y[2][2] + lambda[1]*(y[1][2] - y[2][2]);
                dstsq = CDYN_DIST_DOT(zsol,zsol);
                CDYN_DIST_COPY(y[0],y[2]);
                del[1][0] = del[2][1];
                del[0][0] = del[2][2];
                return(dstsq);
            }
            break;
        }  // END case three points

        case 4: 
        { 
            //
            // case of four points.
            //

            // check optimality of vertex 1.
            d2[2] = del[0][0] - del[1][0];
            d3[4] = del[0][0] - del[2][0];
            d4[8] = del[0][0] - del[3][0];
            if (d2[2]<=0.0f && d3[4]<=0.0f && d4[8]<=0.0f) 
            {
                *nov = 1;
                lambda[0] = d1[0];
                CDYN_DIST_COPY(zsol,y[0]);
                dstsq = del[0][0];
                return(dstsq);
            }

            // check optimality of line segment 1-2
            e132 = del[1][0] - del[2][1];
            e142 = del[1][0] - del[3][1];
            d1[2] = del[1][1] - del[1][0];
            d3[6] = d1[2]*d3[4] + d2[2]*e132;
            d4[11] = d1[2]*d4[8] + d2[2]*e142;
            if (d1[2]>0.0f && d2[2]>0.0f &&
                d3[6]<=0.0f && d4[11]<=0.0f) 
            {
                *nov = 2;
                sum = d1[2] + d2[2];
                lambda[0] = d1[2]/sum;
                lambda[1] = 1.0f - lambda[0];
                zsol[0] = y[1][0] + lambda[0]*(y[0][0] - y[1][0]);
                zsol[1] = y[1][1] + lambda[0]*(y[0][1] - y[1][1]);
                zsol[2] = y[1][2] + lambda[0]*(y[0][2] - y[1][2]);
                dstsq = CDYN_DIST_DOT(zsol,zsol);
                return(dstsq);
            }

            // check optimality of line segment 1-3
            e123 = del[2][0] - del[2][1];
            e143 = del[2][0] - del[3][2];
            d1[4] = del[2][2] - del[2][0];
            d2[6] = d1[4]*d2[2] + d3[4]*e123;
            d4[12] = d1[4]*d4[8] + d3[4]*e143;
            if (d1[4]>0.0f && d2[6]<=0.0f &&
                d3[4]>0.0f && d4[12]<=0.0f) 
            {
                *nov = 2;
                v_index[0][1] = v_index[0][2];
                v_index[1][1] = v_index[1][2];
                sum = d1[4] + d3[4];
                lambda[0] = d1[4]/sum;
                lambda[1] = 1.0f - lambda[0];
                zsol[0] = y[2][0] + lambda[0]*(y[0][0] - y[2][0]);
                zsol[1] = y[2][1] + lambda[0]*(y[0][1] - y[2][1]);
                zsol[2] = y[2][2] + lambda[0]*(y[0][2] - y[2][2]);
                dstsq = CDYN_DIST_DOT(zsol,zsol);
                CDYN_DIST_COPY(y[1],y[2]);
                del[1][0] = del[2][0];
                del[1][1] = del[2][2];
                return(dstsq);
            }

            // check optimality of face 1-2-3
            d2[5] = del[2][2] - del[2][1];
            d3[5] = del[1][1] - del[2][1];
            e213 = -e123;
            d1[6] = d2[5]*d1[2] + d3[5]*e213;
            d4[14] = d1[6]*d4[8] + d2[6]*e142 + d3[6]*e143;
            if (d1[6]>0.0f && d2[6]>0.0f &&
                d3[6]>0.0f && d4[14]<=0.0f) 
            {
                *nov = 3;
                sum = d1[6] + d2[6] + d3[6];
                lambda[0] = d1[6]/sum;
                lambda[1] = d2[6]/sum;
                lambda[2] = 1.0f - lambda[0] - lambda[1];
                zsol[0] = y[2][0] + lambda[0]*(y[0][0] - y[2][0]) + lambda[1]*(y[1][0] - y[2][0]);
                zsol[1] = y[2][1] + lambda[0]*(y[0][1] - y[2][1]) + lambda[1]*(y[1][1] - y[2][1]);
                zsol[2] = y[2][2] + lambda[0]*(y[0][2] - y[2][2]) + lambda[1]*(y[1][2] - y[2][2]);
                dstsq = CDYN_DIST_DOT(zsol,zsol);
                return(dstsq);
            }
            
            // check optimality of line segment 1-4
            e124 = del[3][0] - del[3][1];
            e134 = del[3][0] - del[3][2];
            d1[8] = del[3][3] - del[3][0];
            d2[11] = d1[8]*d2[2] + d4[8]*e124;
            d3[12] = d1[8]*d3[4] + d4[8]*e134;
            if (d1[8]>0.0f && d2[11]<=0.0f &&
                d3[12]<=0.0f && d4[8]>0.0f) 
            {
                *nov = 2;
                v_index[0][1] = v_index[0][3];
                v_index[1][1] = v_index[1][3];
                sum = d1[8] + d4[8];
                lambda[0] = d1[8]/sum;
                lambda[1] = 1.0f - lambda[0];
                zsol[0] = y[3][0] + lambda[0]*(y[0][0] - y[3][0]);
                zsol[1] = y[3][1] + lambda[0]*(y[0][1] - y[3][1]);
                zsol[2] = y[3][2] + lambda[0]*(y[0][2] - y[3][2]);
                dstsq = CDYN_DIST_DOT(zsol,zsol);
                CDYN_DIST_COPY(y[1],y[3]);
                del[1][0] = del[3][0];
                del[1][1] = del[3][3];
                return(dstsq);
            }

            // check optimality of face 1-2-4
            d2[9] = del[3][3] - del[3][1];
            d4[9] = del[1][1] - del[3][1];
            e214 = -e124;
            d1[11] = d2[9]*d1[2] + d4[9]*e214;
            d3[14] = d1[11]*d3[4] + d2[11]*e132 + d4[11]*e134;
            if (d1[11]>0.0f && d2[11]>0.0f &&
                d3[14]<=0.0f && d4[11]>0.0f) 
            {
                *nov = 3;
                v_index[0][2] = v_index[0][3];
                v_index[1][2] = v_index[1][3];
                sum = d1[11] + d2[11] + d4[11];
                lambda[0] = d1[11]/sum;
                lambda[1] = d2[11]/sum;
                lambda[2] = 1.0f - lambda[0] - lambda[1];
                zsol[0] = y[3][0] + lambda[0]*(y[0][0] - y[3][0]) + lambda[1]*(y[1][0] - y[3][0]);
                zsol[1] = y[3][1] + lambda[0]*(y[0][1] - y[3][1]) + lambda[1]*(y[1][1] - y[3][1]);
                zsol[2] = y[3][2] + lambda[0]*(y[0][2] - y[3][2]) + lambda[1]*(y[1][2] - y[3][2]);
                dstsq = CDYN_DIST_DOT(zsol,zsol);
                CDYN_DIST_COPY(y[2],y[3]);
                del[2][0] = del[3][0];
                del[2][1] = del[3][1];
                del[2][2] = del[3][3];
                return(dstsq);
            }

            // check optimality of face 1-3-4
            d3[10] = del[3][3] - del[3][2];
            d4[10] = del[2][2] - del[3][2];
            e314 = -e134;
            d1[12] = d3[10]*d1[4] + d4[10]*e314;
            d2[14] = d1[12]*d2[2] + d3[12]*e123 + d4[12]*e124;
            if (d1[12]>0.0f && d2[14]<=0.0f &&
                d3[12]>0.0f && d4[12]>0.0f) 
            {
                *nov = 3;
                v_index[0][1] = v_index[0][3];
                v_index[1][1] = v_index[1][3];
                sum = d1[12] + d3[12] + d4[12];
                lambda[0] = d1[12]/sum;
                lambda[2] = d3[12]/sum;
                lambda[1] = 1.0f - lambda[0] - lambda[2];
                zsol[0] = y[3][0] + lambda[0]*(y[0][0] - y[3][0]) + lambda[2]*(y[2][0] - y[3][0]);
                zsol[1] = y[3][1] + lambda[0]*(y[0][1] - y[3][1]) + lambda[2]*(y[2][1] - y[3][1]);
                zsol[2] = y[3][2] + lambda[0]*(y[0][2] - y[3][2]) + lambda[2]*(y[2][2] - y[3][2]);
                dstsq = CDYN_DIST_DOT(zsol,zsol);
                CDYN_DIST_COPY(y[1],y[3]);
                del[1][0] = del[3][0];
                del[1][1] = del[3][3];
                del[2][1] = del[3][2];
                return(dstsq);
            }
            
            // check optimality of the hull of all 4 points.
            e243 = del[2][1] - del[3][2];
            d4[13] = d2[5]*d4[9] + d3[5]*e243;
            e234 = del[3][1] - del[3][2];
            d3[13] = d2[9]*d3[5] + d4[9]*e234;
            e324 = -e234;
            d2[13] = d3[10]*d2[5] + d4[10]*e324;
            d1[14] = d2[13]*d1[2] + d3[13]*e213 + d4[13]*e214;
            if (d1[14]>0.0f && d2[14]>0.0f &&
                d3[14]>0.0f && d4[14]>0.0f) 
            {
                sum = d1[14] + d2[14] + d3[14] + d4[14];
                lambda[0] = d1[14]/sum;
                lambda[1] = d2[14]/sum;
                lambda[2] = d3[14]/sum;
                lambda[3] = 1.0f - lambda[0] - lambda[1] - lambda[2];
                zsol[0] = lambda[0]*y[0][0] + lambda[1]*y[1][0] + lambda[2]*y[2][0] + lambda[3]*y[3][0];
                zsol[1] = lambda[0]*y[0][1] + lambda[1]*y[1][1] + lambda[2]*y[2][1] + lambda[3]*y[3][1];
                zsol[2] = lambda[0]*y[0][2] + lambda[1]*y[1][2] + lambda[2]*y[2][2] + lambda[3]*y[3][2];
                dstsq = CDYN_DIST_DOT(zsol,zsol);
                return(dstsq);
            }

            // check optimality of vertex 2.
            if (d1[2]<=0.0f && d3[5]<=0.0f && d4[9]<=0.0f) 
            {
                *nov = 1;
                v_index[0][0] = v_index[0][1];
                v_index[1][0] = v_index[1][1];
                lambda[0] = d2[1];
                CDYN_DIST_COPY(zsol,y[1]);
                dstsq = del[1][1];
                CDYN_DIST_COPY(y[0],y[1]);
                del[0][0] = del[1][1];
                return(dstsq);
            }
            
            // check optimality of vertex 3
            if (d1[4]<=0.0f && d2[5]<=0.0f && d4[10]<=0.0f) 
            {
                *nov = 1;
                v_index[0][0] = v_index[0][2];
                v_index[1][0] = v_index[1][2];
                lambda[0] = d3[3];
                CDYN_DIST_COPY(zsol,y[2]);
                dstsq = del[2][2];
                CDYN_DIST_COPY(y[0],y[2]);
                del[0][0] = del[2][2];
                return(dstsq);
            }
            
            // check optimality of vertex 4
            if (d1[8]<=0.0f && d2[9]<=0.0f && d3[10]<=0.0f) 
            {
                *nov = 1;
                v_index[0][0] = v_index[0][3];
                v_index[1][0] = v_index[1][3];
                lambda[0] = d4[7];
                CDYN_DIST_COPY(zsol,y[3]);
                dstsq = del[3][3];
                CDYN_DIST_COPY(y[0],y[3]);
                del[0][0] = del[3][3];
                return(dstsq);
            }
               
            // check optimality of line segment 2-3
            if (d1[6]<=0.0f && d2[5]>0.0f &&
                d3[5]>0.0f && d4[13]<=0.0f) 
            {
                *nov = 2;
                v_index[0][0] = v_index[0][2];
                v_index[1][0] = v_index[1][2];
                sum = d2[5] + d3[5];
                lambda[1] = d2[5]/sum;
                lambda[0] = 1.0f - lambda[1];
                zsol[0] = y[2][0] + lambda[1]*(y[1][0] - y[2][0]);
                zsol[1] = y[2][1] + lambda[1]*(y[1][1] - y[2][1]);
                zsol[2] = y[2][2] + lambda[1]*(y[1][2] - y[2][2]);
                dstsq = CDYN_DIST_DOT(zsol,zsol);
                CDYN_DIST_COPY(y[0],y[2]);
                del[1][0] = del[2][1];
                del[0][0] = del[2][2];
                return(dstsq);
            }
            
            // check optimality of line segment 2-4
            if (d1[11]<=0.0f && d2[9]>0.0f &&
                d3[13]<=0.0f && d4[9]>0.0f) 
            {
                *nov = 2;
                v_index[0][0] = v_index[0][3];
                v_index[1][0] = v_index[1][3];
                sum = d2[9] + d4[9];
                lambda[1] = d2[9]/sum;
                lambda[0] = 1.0f - lambda[1];
                zsol[0] = y[3][0] + lambda[1]*(y[1][0] - y[3][0]);
                zsol[1] = y[3][1] + lambda[1]*(y[1][1] - y[3][1]);
                zsol[2] = y[3][2] + lambda[1]*(y[1][2] - y[3][2]);
                dstsq = CDYN_DIST_DOT(zsol,zsol);
                CDYN_DIST_COPY(y[0],y[3]);
                del[1][0] = del[3][1];
                del[0][0] = del[3][3];
                return(dstsq);
            }

            // check optimality of line segment 3-4
            if (d1[12]<=0.0f && d2[13]<=0.0f &&
                d3[10]>0.0f && d4[10]>0.0f) 
            {
                *nov = 2;
                v_index[0][0] = v_index[0][2];
                v_index[0][1] = v_index[0][3];
                v_index[1][0] = v_index[1][2];
                v_index[1][1] = v_index[1][3];
                sum = d3[10] + d4[10];
                lambda[0] = d3[10]/sum;
                lambda[1] = 1.0f - lambda[0];
                zsol[0] = y[3][0] + lambda[0]*(y[2][0] - y[3][0]);
                zsol[1] = y[3][1] + lambda[0]*(y[2][1] - y[3][1]);
                zsol[2] = y[3][2] + lambda[0]*(y[2][2] - y[3][2]);
                dstsq = CDYN_DIST_DOT(zsol,zsol);
                CDYN_DIST_COPY(y[0],y[2]);
                CDYN_DIST_COPY(y[1],y[3]);
                del[0][0] = del[2][2];
                del[1][0] = del[3][2];
                del[1][1] = del[3][3];
                return(dstsq);
            }
             
            // check optimality of face 2-3-4
            if (d1[14]<=0.0f && d2[13]>0.0f &&
                d3[13]>0.0f && d4[13]>0.0f) 
            {
                *nov = 3;
                v_index[0][0] = v_index[0][3];
                v_index[1][0] = v_index[1][3];
                sum = d2[13] + d3[13] + d4[13];
                lambda[1] = d2[13]/sum;
                lambda[2] = d3[13]/sum;
                lambda[0] = 1.0f - lambda[1] - lambda[2];
                zsol[0] = y[3][0] + lambda[1]*(y[1][0] - y[3][0]) + lambda[2]*(y[2][0] - y[3][0]);
                zsol[1] = y[3][1] + lambda[1]*(y[1][1] - y[3][1]) + lambda[2]*(y[2][1] - y[3][1]);
                zsol[2] = y[3][2] + lambda[1]*(y[1][2] - y[3][2]) + lambda[2]*(y[2][2] - y[3][2]);
                dstsq = CDYN_DIST_DOT(zsol,zsol);
                CDYN_DIST_COPY(y[0],y[3]);
                del[0][0] = del[3][3];
                del[1][0] = del[3][1];
                del[2][0] = del[3][2];
                return(dstsq);
            }
            break;
        } // END case four points

        default: 
        {
            // cDynPrintf(stderr,"Invalid value for nov %d given \n",*nov);
            break;
        }
    } /* END switch */
    }

  //----------------------------------------------------------------------
  //  The  backup procedure  begins ...                                  
  //----------------------------------------------------------------------

  switch (*nov) 
  {
    case 1: 
    { 
        //
        // case of a single point.
        //
        dstsq = del[0][0];
        lambda[0] = d1[0];
        CDYN_DIST_COPY(zsol,y[0]);
        return(dstsq);
    }  // END case of a single point

    case 2: 
    { 
        //
        // case of a two points.
        // 

        if (backup) 
        {
            d2[2] = del[0][0] - del[1][0];
            d1[2] = del[1][1] - del[1][0];
        }

        // check vertex 1
        dstsq = del[0][0];
        novd = 1;
        lambda[0] = d1[0];
        CDYN_DIST_COPY(zsol,y[0]);
        iord[0] = 0;
      
        // check line segment 1-2
        if (d1[2]>0.0f && d2[2]>0.0f) 
        {
            sum = d1[2] + d2[2];
            alsd[0] = d1[2]/sum;
            alsd[1] = 1.0f - alsd[0];
            zsold[0] = y[1][0] + alsd[0]*(y[0][0] - y[1][0]);
            zsold[1] = y[1][1] + alsd[0]*(y[0][1] - y[1][1]);
            zsold[2] = y[1][2] + alsd[0]*(y[0][2] - y[1][2]);
            dstsqd = CDYN_DIST_DOT(zsold,zsold);
            if (dstsqd < dstsq) 
            {
                dstsq = dstsqd;
                novd = 2;
                lambda[0] = alsd[0];
                lambda[1] = alsd[1];
                CDYN_DIST_COPY(zsol,zsold);
                iord[0] = 0;
                iord[1] = 1;
            }
        }

        // check vertex 2
        if (del[1][1] < dstsq) 
        {
            dstsq = del[1][1];
            novd = 1;
            lambda[0] = d2[1];
            CDYN_DIST_COPY(zsol,y[1]);
            iord[0] = 1;
        }
        break;
    }  // END case two points.


    case 3: 
    { 
        //
        // case of a three points.
        // 
     
        if (backup) 
        {
            d2[2] = del[0][0] - del[1][0];
            d3[4] = del[0][0] - del[2][0];
            e132 = del[1][0] - del[2][1];
            d1[2] = del[1][1] - del[1][0];
            d3[6] = d1[2]*d3[4] + d2[2]*e132;
            e123 = del[2][0] - del[2][1];
            d1[4] = del[2][2] - del[2][0];
            d2[6] = d1[4]*d2[2] + d3[4]*e123;
            e213 = -e123;
            d2[5] = del[2][2] - del[2][1];
            d3[5] = del[1][1] - del[2][1];
            d1[6] = d2[5]*d1[2] + d3[5]*e213;
        }

        // check vertex 1
        dstsq = del[0][0];
        novd = 1;
        lambda[0] = d1[0];
        CDYN_DIST_COPY(zsol,y[0]);
        iord[0] = 0;

        // check line segment 1-2
        if (d1[2]>0.0f && d2[2]>0.0f) 
        {
            sum = d1[2] + d2[2];
            alsd[0] = d1[2]/sum;
            alsd[1] = 1.0f - alsd[0];
            zsold[0] = y[1][0] + alsd[0]*(y[0][0] - y[1][0]);
            zsold[1] = y[1][1] + alsd[0]*(y[0][1] - y[1][1]);
            zsold[2] = y[1][2] + alsd[0]*(y[0][2] - y[1][2]);
            dstsqd = CDYN_DIST_DOT(zsold,zsold);
            if (dstsqd < dstsq) 
            {
                dstsq = dstsqd;
                novd = 2;
                lambda[0] = alsd[0];
                lambda[1] = alsd[1];
                CDYN_DIST_COPY(zsol,zsold);
                iord[0] = 0;
                iord[1] = 1;
            }
        }

        // check line segment 1-3
        if (d1[4]>0.0f && d3[4]>0.0f) 
        {
            sum = d1[4] + d3[4];
            alsd[0] = d1[4]/sum;
            alsd[1] = 1.0f - alsd[0];
            zsold[0] = y[2][0] + alsd[0]*(y[0][0] - y[2][0]);
            zsold[1] = y[2][1] + alsd[0]*(y[0][1] - y[2][1]);
            zsold[2] = y[2][2] + alsd[0]*(y[0][2] - y[2][2]);
            dstsqd = CDYN_DIST_DOT(zsold,zsold);
            if (dstsqd < dstsq) 
            {
                dstsq = dstsqd;
                novd = 2;
                lambda[0] = alsd[0];
                lambda[1] = alsd[1];
                CDYN_DIST_COPY(zsol,zsold);
                iord[0] = 0;
                iord[1] = 2;
            }
        }

        // check face 1-2-3
        if (d1[6]>0.0f && d2[6]>0.0f && d3[6]>0.0f) 
        {
            sum = d1[6] + d2[6] + d3[6];
            alsd[0] = d1[6]/sum;
            alsd[1] = d2[6]/sum;
            alsd[2] = 1.0f - alsd[0] - alsd[1];
            zsold[0] = y[2][0] + alsd[0]*(y[0][0] - y[2][0]) + alsd[1]*(y[1][0] - y[2][0]);
            zsold[1] = y[2][1] + alsd[0]*(y[0][1] - y[2][1]) + alsd[1]*(y[1][1] - y[2][1]);
            zsold[2] = y[2][2] + alsd[0]*(y[0][2] - y[2][2]) + alsd[1]*(y[1][2] - y[2][2]);
            dstsqd = CDYN_DIST_DOT(zsold,zsold);
            if (dstsqd < dstsq) 
            {
                dstsq = dstsqd;
                novd = 3;
                lambda[0] = alsd[0];
                lambda[1] = alsd[1];
                lambda[2] = alsd[2];
                CDYN_DIST_COPY(zsol,zsold);
                iord[0] = 0;
                iord[1] = 1;
                iord[2] = 2;
            }
        }

        // check vertex 2
        if (del[1][1] < dstsq) 
        {
            novd = 1;
            dstsq = del[1][1];
            lambda[0] = d2[1];
            CDYN_DIST_COPY(zsol,y[1]);
            iord[0] = 1;
        }

        // check vertex 3
        if (del[2][2] < dstsq) 
        {
            novd = 1;
            dstsq = del[2][2];
            lambda[0] = d3[3];
            CDYN_DIST_COPY(zsol,y[2]);
            iord[0] = 2;
        }

        // check line segment 2-3
        if (d2[5]>0.0f && d3[5]>0.0f) 
        {
            sum = d2[5] + d3[5];
            alsd[1] = d2[5]/sum;
            alsd[0] = 1.0f - alsd[1];
            zsold[0] = y[2][0] + alsd[1]*(y[1][0] - y[2][0]);
            zsold[1] = y[2][1] + alsd[1]*(y[1][1] - y[2][1]);
            zsold[2] = y[2][2] + alsd[1]*(y[1][2] - y[2][2]);
            dstsqd = CDYN_DIST_DOT(zsold,zsold);
            if (dstsqd < dstsq) 
            {
                dstsq = dstsqd;
                novd = 2;
                lambda[0] = alsd[0];
                lambda[1] = alsd[1];
                CDYN_DIST_COPY(zsol,zsold);
                iord[0] = 2;
                iord[1] = 1;
            }
        }
        break;
        }  // END case three points


    case 4: 
    { 
        //
        // case of a four points.
        // 

        if (backup) 
        {
            d2[2] = del[0][0] - del[1][0];
            d3[4] = del[0][0] - del[2][0];
            d4[8] = del[0][0] - del[3][0];
            e132 = del[1][0] - del[2][1];
            e142 = del[1][0] - del[3][1];
            d1[2] = del[1][1] - del[1][0];
            d3[6] = d1[2]*d3[4] + d2[2]*e132;
            d4[11] = d1[2]*d4[8] + d2[2]*e142;
            e123 = del[2][0] - del[2][1];
            e143 = del[2][0] - del[3][2];
            d1[4] = del[2][2] - del[2][0];
            d2[6] = d1[4]*d2[2] + d3[4]*e123;
            d4[12] = d1[4]*d4[8] + d3[4]*e143;
            d2[5] = del[2][2] - del[2][1];
            d3[5] = del[1][1] - del[2][1];
            e213 = -e123;
            d1[6] = d2[5]*d1[2] + d3[5]*e213;
            d4[14] = d1[6]*d4[8] + d2[6]*e142 + d3[6]*e143;
            e124 = del[3][0] - del[3][1];
            e134 = del[3][0] - del[3][2];
            d1[8] = del[3][3] - del[3][0];
            d2[11] = d1[8]*d2[2] + d4[8]*e124;
            d3[12] = d1[8]*d3[4] + d4[8]*e134;
            d2[9] = del[3][3] - del[3][1];
            d4[9] = del[1][1] - del[3][1];
            e214 = -e124;
            d1[11] = d2[9]*d1[2] + d4[9]*e214;
            d3[14] = d1[11]*d3[4] + d2[11]*e132 + d4[11]*e134;
            d3[10] = del[3][3] - del[3][2];
            d4[10] = del[2][2] - del[3][2];
            e314 = -e134;
            d1[12] = d3[10]*d1[4] + d4[10]*e314;
            d2[14] = d1[12]*d2[2] + d3[12]*e123 + d4[12]*e124;
            e243 = del[2][1] - del[3][2];
            d4[13] = d2[5]*d4[9] + d3[5]*e243;
            e234 = del[3][1] - del[3][2];
            d3[13] = d2[9]*d3[5] + d4[9]*e234;
            e324 = -e234;
            d2[13] = d3[10]*d2[5] + d4[10]*e324;
            d1[14] = d2[13]*d1[2] + d3[13]*e213 + d4[13]*e214;
        }

        // check vertex 1
        dstsq = del[0][0];
        novd = 1;
        lambda[0] = d1[0];
        CDYN_DIST_COPY(zsol,y[0]);
        iord[0] = 0;

        // check line segment 1-2
        if (d1[2]>0.0f && d2[2]>0.0f) 
        {
            sum = d1[2] + d2[2];
            alsd[0] = d1[2]/sum;
            alsd[1] = 1.0f - alsd[0];
            zsold[0] = y[1][0] + alsd[0]*(y[0][0] - y[1][0]);
            zsold[1] = y[1][1] + alsd[0]*(y[0][1] - y[1][1]);
            zsold[2] = y[1][2] + alsd[0]*(y[0][2] - y[1][2]);
            dstsqd = CDYN_DIST_DOT(zsold,zsold);
            if (dstsqd < dstsq) {
                dstsq = dstsqd;
                novd = 2;
                lambda[0] = alsd[0];
                lambda[1] = alsd[1];
                CDYN_DIST_COPY(zsol,zsold);
                iord[0] = 0;
                iord[1] = 1;
            }
        }

        // check line segment 1-3
        if (d1[4]>0.0f && d3[4]>0.0f) 
        {
            sum = d1[4] + d3[4];
            alsd[0] = d1[4]/sum;
            alsd[1] = 1.0f - alsd[0];
            zsold[0] = y[2][0] + alsd[0]*(y[0][0] - y[2][0]);
            zsold[1] = y[2][1] + alsd[0]*(y[0][1] - y[2][1]);
            zsold[2] = y[2][2] + alsd[0]*(y[0][2] - y[2][2]);
            dstsqd = CDYN_DIST_DOT(zsold,zsold);
            if (dstsqd < dstsq) 
            {
                dstsq = dstsqd;
                novd = 2;
                lambda[0] = alsd[0];
                lambda[1] = alsd[1];
                CDYN_DIST_COPY(zsol,zsold);
                iord[0] = 0;
                iord[1] = 2;
            }
        }

        // check face 1-2-3
        if (d1[6]>0.0f && d2[6]>0.0f && d3[6]>0.0f) 
        {
            sum = d1[6] + d2[6] + d3[6];
            alsd[0] = d1[6]/sum;
            alsd[1] = d2[6]/sum;
            alsd[2] = 1.0f - alsd[0] - alsd[1];
            zsold[0] = y[2][0] + alsd[0]*(y[0][0] - y[2][0]) + alsd[1]*(y[1][0] - y[2][0]);
            zsold[1] = y[2][1] + alsd[0]*(y[0][1] - y[2][1]) + alsd[1]*(y[1][1] - y[2][1]);
            zsold[2] = y[2][2] + alsd[0]*(y[0][2] - y[2][2]) + alsd[1]*(y[1][2] - y[2][2]);
            dstsqd = CDYN_DIST_DOT(zsold,zsold);
            if (dstsqd < dstsq) 
            {
                dstsq = dstsqd;
                novd = 3;
                lambda[0] = alsd[0];
                lambda[1] = alsd[1];
                lambda[2] = alsd[2];
                CDYN_DIST_COPY(zsol,zsold);
                iord[0] = 0;
                iord[1] = 1;
                iord[2] = 2;
            }
        }

      // check line segment 1-4
        if (d1[8]>0.0f && d4[8]>0.0f) 
        {
            sum = d1[8] + d4[8];
            alsd[0] = d1[8]/sum;
            alsd[1] = 1.0f - alsd[0];
            zsold[0] = y[3][0] + alsd[0]*(y[0][0] - y[3][0]);
            zsold[1] = y[3][1] + alsd[0]*(y[0][1] - y[3][1]);
            zsold[2] = y[3][2] + alsd[0]*(y[0][2] - y[3][2]);
            dstsqd = CDYN_DIST_DOT(zsold,zsold);
            if (dstsqd < dstsq) 
            {
                dstsq = dstsqd;
                novd = 2;
                lambda[0] = alsd[0];
                lambda[1] = alsd[1];
                CDYN_DIST_COPY(zsol,zsold);
                iord[0] = 0;
                iord[1] = 3;
            }
        }
      
        // check face 1-2-4
        if (d1[11]>0.0f && d2[11]>0.0f && d4[11]>0.0f) 
        {
            sum = d1[11] + d2[11] + d4[11];
            alsd[0] = d1[11]/sum;
            alsd[1] = d2[11]/sum;
            alsd[2] = 1.0f - alsd[0] - alsd[1];
            zsold[0] = y[3][0] + alsd[0]*(y[0][0] - y[3][0]) + alsd[1]*(y[1][0] - y[3][0]);
            zsold[1] = y[3][1] + alsd[0]*(y[0][1] - y[3][1]) + alsd[1]*(y[1][1] - y[3][1]);
            zsold[2] = y[3][2] + alsd[0]*(y[0][2] - y[3][2]) + alsd[1]*(y[1][2] - y[3][2]);
            dstsqd = CDYN_DIST_DOT(zsold,zsold);
            if (dstsqd < dstsq) 
            {
                dstsq = dstsqd;
                novd = 3;
                lambda[0] = alsd[0];
                lambda[1] = alsd[1];
                lambda[2] = alsd[2];
                CDYN_DIST_COPY(zsol,zsold);
                iord[0] = 0;
                iord[1] = 1;
                iord[2] = 3;
            }
        }

       // check face 1-3-4
        if (d1[12]>0.0f && d3[12]>0.0f && d4[12]>0.0f) 
        {
            sum = d1[12] + d3[12] + d4[12];
            alsd[0] = d1[12]/sum;
            alsd[2] = d3[12]/sum;
            alsd[1] = 1.0f - alsd[0] - alsd[2];
            zsold[0] = y[3][0] + alsd[0]*(y[0][0] - y[3][0]) + alsd[2]*(y[2][0] - y[3][0]);
            zsold[1] = y[3][1] + alsd[0]*(y[0][1] - y[3][1]) + alsd[2]*(y[2][1] - y[3][1]);
            zsold[2] = y[3][2] + alsd[0]*(y[0][2] - y[3][2]) + alsd[2]*(y[2][2] - y[3][2]);
            dstsqd = CDYN_DIST_DOT(zsold,zsold);
            if (dstsqd < dstsq) {
                dstsq = dstsqd;
                novd = 3;
                lambda[0] = alsd[0];
                lambda[1] = alsd[1];
                lambda[2] = alsd[2];
                CDYN_DIST_COPY(zsol,zsold);
                iord[0] = 0;
                iord[1] = 3;
                iord[2] = 2;
            }
        }

        // check the hull of all 4 points
        if (d1[14]>0.0f && d2[14]>0.0f &&
            d3[14]>0.0f && d4[14]>0.0f) 
        {
                sum = d1[14] + d2[14] + d3[14] + d4[14];
                alsd[0] = d1[14]/sum;
                alsd[1] = d2[14]/sum;
                alsd[2] = d3[14]/sum;
                alsd[3] = 1.0f - alsd[0] - alsd[1] - alsd[2];
                zsold[0] = alsd[0]*y[0][0] + alsd[1]*y[1][0] + alsd[2]*y[2][0] + alsd[3]*y[3][0];
                zsold[1] = alsd[0]*y[0][1] + alsd[1]*y[1][1] + alsd[2]*y[2][1] + alsd[3]*y[3][1];
                zsold[2] = alsd[0]*y[0][2] + alsd[1]*y[1][2] + alsd[2]*y[2][2] + alsd[3]*y[3][2];
                dstsqd = CDYN_DIST_DOT(zsold,zsold);
                if (dstsqd < dstsq) 
                {
                    dstsq = dstsqd;
                    novd = 4;
                    lambda[0] = alsd[0];
                    lambda[1] = alsd[1];
                    lambda[2] = alsd[2];
                    lambda[3] = alsd[3];
                    CDYN_DIST_COPY(zsol,zsold);
                    iord[0] = 0;
                    iord[1] = 1;
                    iord[2] = 2;
                    iord[3] = 3;
                }
        }
      
        // check vertex 2
        if (del[1][1] < dstsq) 
        {
            novd = 1;
            dstsq = del[1][1];
            lambda[0] = d2[1];
            CDYN_DIST_COPY(zsol,y[1]);
            iord[0] = 1;
        }

        // check vertex 3
        if (del[2][2] < dstsq) 
        {
            novd = 1;
            dstsq = del[2][2];
            lambda[0] = d3[3];
            CDYN_DIST_COPY(zsol,y[2]);
            iord[0] = 2;
        }

        // check vertex 4
        if (del[3][3] < dstsq) 
        {
            novd = 1;
            dstsq = del[3][3];
            lambda[0] = d4[7];
            CDYN_DIST_COPY(zsol,y[3]);
            iord[0] = 3;
        }

        // check line segment 2-3
        if (d2[5]>0.0f && d3[5]>0.0f) 
        {
            sum = d2[5] + d3[5];
            alsd[1] = d2[5]/sum;
            alsd[0] = 1.0f - alsd[1];
            zsold[0] = y[2][0] + alsd[1]*(y[1][0] - y[2][0]);
            zsold[1] = y[2][1] + alsd[1]*(y[1][1] - y[2][1]);
            zsold[2] = y[2][2] + alsd[1]*(y[1][2] - y[2][2]);
            dstsqd = CDYN_DIST_DOT(zsold,zsold);
            if (dstsqd < dstsq) {
                dstsq = dstsqd;
                novd = 2;
                lambda[0] = alsd[0];
                lambda[1] = alsd[1];
                CDYN_DIST_COPY(zsol,zsold);
                iord[0] = 2;
                iord[1] = 1;
            }
        }

        // check line segment 2-4
        if (d2[9]>0.0f && d4[9]>0.0f) 
        {
            sum = d2[9] + d4[9];
            alsd[1] = d2[9]/sum;
            alsd[0] = 1.0f - alsd[1];
            zsold[0] = y[3][0] + alsd[1]*(y[1][0] - y[3][0]);
            zsold[1] = y[3][1] + alsd[1]*(y[1][1] - y[3][1]);
            zsold[2] = y[3][2] + alsd[1]*(y[1][2] - y[3][2]);
            dstsqd = CDYN_DIST_DOT(zsold,zsold);
            if (dstsqd < dstsq) 
            {
                dstsq = dstsqd;
                novd = 2;
                lambda[0] = alsd[0];
                lambda[1] = alsd[1];
                CDYN_DIST_COPY(zsol,zsold);
                iord[0] = 3;
                iord[1] = 1;
            }
        }

        // check line segment 3-4
        if (d3[10]>0.0f && d4[10]>0.0f) 
        {
            sum = d3[10] + d4[10];
            alsd[0] = d3[10]/sum;
            alsd[1] = 1.0f - alsd[0];
            zsold[0] = y[3][0] + alsd[0]*(y[2][0] - y[3][0]);
            zsold[1] = y[3][1] + alsd[0]*(y[2][1] - y[3][1]);
            zsold[2] = y[3][2] + alsd[0]*(y[2][2] - y[3][2]);
            dstsqd = CDYN_DIST_DOT(zsold,zsold);
            if (dstsqd < dstsq) {
                dstsq = dstsqd;
                novd = 2;
                lambda[0] = alsd[0];
                lambda[1] = alsd[1];
                CDYN_DIST_COPY(zsol,zsold);
                iord[0] = 2;
                iord[1] = 3;
            }
        }

        /* check face 2-3-4 */
        if (d2[13]>0.0f && d3[13]>0.0f && d4[13]>0.0f) 
        {
            sum = d2[13] + d3[13] + d4[13];
            alsd[1] = d2[13]/sum;
            alsd[2] = d3[13]/sum;
            alsd[0] = 1.0f - alsd[1] - alsd[2];
            zsold[0] = y[3][0] + alsd[1]*(y[1][0] - y[3][0]) + alsd[2]*(y[2][0] - y[3][0]);
            zsold[1] = y[3][1] + alsd[1]*(y[1][1] - y[3][1]) + alsd[2]*(y[2][1] - y[3][1]);
            zsold[2] = y[3][2] + alsd[1]*(y[1][2] - y[3][2]) + alsd[2]*(y[2][2] - y[3][2]);
            dstsqd = CDYN_DIST_DOT(zsold,zsold);
            if (dstsqd < dstsq) 
            {
                dstsq = dstsqd;
                novd = 3;
                lambda[0] = alsd[0];
                lambda[1] = alsd[1];
                lambda[2] = alsd[2];
                CDYN_DIST_COPY(zsol,zsold);
                iord[0] = 3;
                iord[1] = 1;
                iord[2] = 2;
            }
        }
        break;
        
    } // END case of four points.
    } // END switch

    // final reordering
    for (k=0;k<*nov;k++) 
    {
        v_indexd[0][k]=v_index[0][k];
        v_indexd[1][k]=v_index[1][k];
        CDYN_DIST_COPY(yd[k],y[k]);

        for (l=0;l<k;l++) 
        {
            deld[k][l]=del[k][l];
        }
        deld[k][k]=del[k][k];
    }

    *nov = novd;
    for (k=0;k<*nov;k++) 
    {
        kk = iord[k];
        v_index[0][k] = v_indexd[0][kk];
        v_index[1][k] = v_indexd[1][kk];
        CDYN_DIST_COPY(y[k],yd[kk]);
        for (l=0;l<k;l++) 
        {
            ll = iord[l];
            if (kk >= ll)
                del[k][l] = deld[kk][ll];
            else
                del[k][l] = deld[ll][kk];
        }
        del[k][k] = deld[kk][kk];
    }

    return(dstsq);
}

//---------------------------------------------------------------------------

static double _cDynDistDSBP2(int *nov,int v_index[2][4],
    double lambda[4],
    double y[4][3],
    double del[4][4],
    double zsol[3])
{
    double sum;
    double e132,e142,e123,e143,e213,e243;
    double e124,e134,e214,e234,e314,e324;
    double d1[15],d2[15],d3[15],d4[15];
    double dstsq=0.0f;

    d1[0]=d2[1]=d3[3]=d4[7]=1.0f;

    //
    // regular distance subalgoritm begins
    //
    switch (*nov) 
    {
        case 1: 
        {       
            //
            // case of a single point.
            //
            lambda[0] = d1[0];
            CDYN_DIST_COPY(zsol,y[0]);
            dstsq = del[0][0];
            return(dstsq);
        } // END case one point

        case 2: 
        {       
            //
            // case of two points.
            //

            // check optimality of vertex 1.
            d2[2] = del[0][0] - del[1][0];
            if (d2[2] <= 0.0f) 
            {
                *nov = 1;
                lambda[0] = 1.0;
                CDYN_DIST_COPY(zsol,y[0]);
                dstsq = del[0][0];
                return(dstsq);
            }
            
            // check optimality of line segment 1-2
            d1[2] = del[1][1] - del[1][0];
            if (d1[2]>0.0f && d2[2]>0.0f) 
            {
                sum = d1[2] + d2[2];
                lambda[0] = d1[2]/sum;
                lambda[1] = 1.0f - lambda[0];
                zsol[0] = y[1][0] + lambda[0]*(y[0][0] - y[1][0]);
                zsol[1] = y[1][1] + lambda[0]*(y[0][1] - y[1][1]);
                zsol[2] = y[1][2] + lambda[0]*(y[0][2] - y[1][2]);
                dstsq = CDYN_DIST_DOT(zsol,zsol);
                return(dstsq);
            }

            // check optimality of vertex 2
            if (d1[2]<=0.0f) 
            {
                *nov = 1;
                v_index[0][0] = v_index[0][1];
                v_index[1][0] = v_index[1][1];
                lambda[0] = d2[1];
                CDYN_DIST_COPY(zsol,y[1]);
                dstsq = del[1][1];
                CDYN_DIST_COPY(y[0],y[1]);
                del[0][0] = del[1][1];
                return(dstsq);
            }
            break;
        }  // END case two points

        
        case 3: 
        {       
            //
            // case of three points.
            //

            // check optimality of vertex 1 */
            d2[2] = del[0][0] - del[1][0];
            d3[4] = del[0][0] - del[2][0];
            if (d2[2]<=0.0f && d3[4]<=0.0f)
            {
                *nov = 1;
                lambda[0] = d1[0];
                CDYN_DIST_COPY(zsol,y[0]);
                dstsq = del[0][0];
                return(dstsq);
            }
                
            // check optimality of line segment 1-2
            e132 = del[1][0] - del[2][1];
            d1[2] = del[1][1] - del[1][0];
            d3[6] = d1[2]*d3[4] + d2[2]*e132;
            if (d1[2]>0.0f && d2[2]> 0.0f && d3[6]<=0.0f) 
            {
                *nov = 2;
                sum = d1[2] + d2[2];
                lambda[0] = d1[2]/sum;
                lambda[1] = 1.0f - lambda[0];
                zsol[0] = y[1][0] + lambda[0]*(y[0][0] - y[1][0]);
                zsol[1] = y[1][1] + lambda[0]*(y[0][1] - y[1][1]);
                zsol[2] = y[1][2] + lambda[0]*(y[0][2] - y[1][2]);
                dstsq = CDYN_DIST_DOT(zsol,zsol);
                return(dstsq);
            }
                
            // check optimality of line segment 1-3
            e123 = del[2][0] - del[2][1];
            d1[4] = del[2][2] - del[2][0];
            d2[6] = d1[4]*d2[2] + d3[4]*e123;
            if (d1[4]>0.0f && d2[6]<=0.0f && d3[4]>0.0f) 
            {
                *nov = 2;
                v_index[0][1] = v_index[0][2];
                v_index[1][1] = v_index[1][2];
                sum = d1[4] + d3[4];
                lambda[0] = d1[4]/sum;
                lambda[1] = 1.0f - lambda[0];
                zsol[0] = y[2][0] + lambda[0]*(y[0][0] - y[2][0]);
                zsol[1] = y[2][1] + lambda[0]*(y[0][1] - y[2][1]);
                zsol[2] = y[2][2] + lambda[0]*(y[0][2] - y[2][2]);
                dstsq = CDYN_DIST_DOT(zsol,zsol);
                CDYN_DIST_COPY(y[1],y[2]);
                del[1][0] = del[2][0];
                del[1][1] = del[2][2];
                return(dstsq);
            }
        
            // check optimality of face 1-2-3
            e213 = -e123;
            d2[5] = del[2][2] - del[2][1];
            d3[5] = del[1][1] - del[2][1];
            d1[6] = d2[5]*d1[2] + d3[5]*e213;
            if (d1[6]>0.0f && d2[6]>0.0f && d3[6]>0.0f) 
            {
                sum = d1[6] + d2[6] + d3[6];
                lambda[0] = d1[6]/sum;
                lambda[1] = d2[6]/sum;
                lambda[2] = 1.0f - lambda[0] - lambda[1];
                zsol[0] = y[2][0] + lambda[0]*(y[0][0] - y[2][0]) + lambda[1]*(y[1][0] - y[2][0]);
                zsol[1] = y[2][1] + lambda[0]*(y[0][1] - y[2][1]) + lambda[1]*(y[1][1] - y[2][1]);
                zsol[2] = y[2][2] + lambda[0]*(y[0][2] - y[2][2]) + lambda[1]*(y[1][2] - y[2][2]);
                dstsq = CDYN_DIST_DOT(zsol,zsol);
                //	  cDynPrintf("case 3: sum=%13.9f\n",sum);
                return(dstsq);
            }
        
            // check optimality of vertex 2
            if (d1[2]<=0.0f && d3[5]<=0.0f) 
            {
                *nov = 1;
                v_index[0][0] = v_index[0][1];
                v_index[1][0] = v_index[1][1];
                lambda[0] = d2[1];
                CDYN_DIST_COPY(zsol,y[1]);
                dstsq = del[1][1];
                CDYN_DIST_COPY(y[0],y[1]);
                del[0][0] = del[1][1];
                return(dstsq);
            }
            
            // check optimality of vertex 3
            if (d1[4]<=0.0f && d2[5]<=0.0f) 
            {
                *nov = 1;
                v_index[0][0] = v_index[0][2];
                v_index[1][0] = v_index[1][2];
                lambda[0] = d3[3];
                CDYN_DIST_COPY(zsol,y[2]);
                dstsq = del[2][2];
                CDYN_DIST_COPY(y[0],y[2]);
                del[0][0] = del[2][2];
                return(dstsq);
            }

            // check optimality of line segment 2-3
            if (d1[6]<=0.0f && d2[5]>0.0f && d3[5]>0.0f) 
            {
                *nov = 2;
                v_index[0][0] = v_index[0][2];
                v_index[1][0] = v_index[1][2];
                sum = d2[5] + d3[5];
                lambda[1] = d2[5]/sum;
                lambda[0] = 1.0f - lambda[1];
                zsol[0] = y[2][0] + lambda[1]*(y[1][0] - y[2][0]);
                zsol[1] = y[2][1] + lambda[1]*(y[1][1] - y[2][1]);
                zsol[2] = y[2][2] + lambda[1]*(y[1][2] - y[2][2]);
                dstsq = CDYN_DIST_DOT(zsol,zsol);
                CDYN_DIST_COPY(y[0],y[2]);
                del[1][0] = del[2][1];
                del[0][0] = del[2][2];
                return(dstsq);
            }
            
            break;
        } // END case three points
        
        case 4: 
        { 
            //
            // case of three points.
            //
            
            // check optimality of vertex 1
            d2[2] = del[0][0] - del[1][0];
            d3[4] = del[0][0] - del[2][0];
            d4[8] = del[0][0] - del[3][0];
            if (d2[2]<=0.0f && d3[4]<=0.0f && d4[8]<=0.0f) 
            {
                *nov = 1;
                lambda[0] = d1[0];
                CDYN_DIST_COPY(zsol,y[0]);
                dstsq = del[0][0];
                return(dstsq);
            }
            
            // check optimality of line segment 1-2
            e132 = del[1][0] - del[2][1];
            e142 = del[1][0] - del[3][1];
            d1[2] = del[1][1] - del[1][0];
            d3[6] = d1[2]*d3[4] + d2[2]*e132;
            d4[11] = d1[2]*d4[8] + d2[2]*e142;
            if (d1[2]>0.0f && d2[2]>0.0f &&
            d3[6]<=0.0f && d4[11]<=0.0f) 
            {
                *nov = 2;
                sum = d1[2] + d2[2];
                lambda[0] = d1[2]/sum;
                lambda[1] = 1.0f - lambda[0];
                zsol[0] = y[1][0] + lambda[0]*(y[0][0] - y[1][0]);
                zsol[1] = y[1][1] + lambda[0]*(y[0][1] - y[1][1]);
                zsol[2] = y[1][2] + lambda[0]*(y[0][2] - y[1][2]);
                dstsq = CDYN_DIST_DOT(zsol,zsol);
                return(dstsq);
            }
        
            // check optimality of line segment 1-3
            e123 = del[2][0] - del[2][1];
            e143 = del[2][0] - del[3][2];
            d1[4] = del[2][2] - del[2][0];
            d2[6] = d1[4]*d2[2] + d3[4]*e123;
            d4[12] = d1[4]*d4[8] + d3[4]*e143;
            if (d1[4]>0.0f && d2[6]<=0.0f &&
                d3[4]>0.0f && d4[12]<=0.0f) 
            {
                    *nov = 2;
                    v_index[0][1] = v_index[0][2];
                    v_index[1][1] = v_index[1][2];
                    sum = d1[4] + d3[4];
                    lambda[0] = d1[4]/sum;
                    lambda[1] = 1.0f - lambda[0];
                    zsol[0] = y[2][0] + lambda[0]*(y[0][0] - y[2][0]);
                    zsol[1] = y[2][1] + lambda[0]*(y[0][1] - y[2][1]);
                    zsol[2] = y[2][2] + lambda[0]*(y[0][2] - y[2][2]);
                    dstsq = CDYN_DIST_DOT(zsol,zsol);
                    CDYN_DIST_COPY(y[1],y[2]);
                    del[1][0] = del[2][0];
                    del[1][1] = del[2][2];
                    return(dstsq);
            }
            
            // check optimality of face 1-2-3
            d2[5] = del[2][2] - del[2][1];
            d3[5] = del[1][1] - del[2][1];
            e213 = -e123;
            d1[6] = d2[5]*d1[2] + d3[5]*e213;
            d4[14] = d1[6]*d4[8] + d2[6]*e142 + d3[6]*e143;
            if (d1[6]>0.0f && d2[6]>0.0f &&
                d3[6]>0.0f && d4[14]<=0.0f) 
            {
                *nov = 3;
                sum = d1[6] + d2[6] + d3[6];
                lambda[0] = d1[6]/sum;
                lambda[1] = d2[6]/sum;
                lambda[2] = 1.0f - lambda[0] - lambda[1];
                zsol[0] = y[2][0] + lambda[0]*(y[0][0] - y[2][0]) + lambda[1]*(y[1][0] - y[2][0]);
                zsol[1] = y[2][1] + lambda[0]*(y[0][1] - y[2][1]) + lambda[1]*(y[1][1] - y[2][1]);
                zsol[2] = y[2][2] + lambda[0]*(y[0][2] - y[2][2]) + lambda[1]*(y[1][2] - y[2][2]);
                dstsq = CDYN_DIST_DOT(zsol,zsol);
                return(dstsq);
            }
            
            // check optimality of line segment 1-4
            e124 = del[3][0] - del[3][1];
            e134 = del[3][0] - del[3][2];
            d1[8] = del[3][3] - del[3][0];
            d2[11] = d1[8]*d2[2] + d4[8]*e124;
            d3[12] = d1[8]*d3[4] + d4[8]*e134;
            if (d1[8]>0.0f && d2[11]<=0.0f &&
                d3[12]<=0.0f && d4[8]>0.0f) 
            {
                *nov = 2;
                v_index[0][1] = v_index[0][3];
                v_index[1][1] = v_index[1][3];
                sum = d1[8] + d4[8];
                lambda[0] = d1[8]/sum;
                lambda[1] = 1.0f - lambda[0];
                zsol[0] = y[3][0] + lambda[0]*(y[0][0] - y[3][0]);
                zsol[1] = y[3][1] + lambda[0]*(y[0][1] - y[3][1]);
                zsol[2] = y[3][2] + lambda[0]*(y[0][2] - y[3][2]);
                dstsq = CDYN_DIST_DOT(zsol,zsol);
                CDYN_DIST_COPY(y[1],y[3]);
                del[1][0] = del[3][0];
                del[1][1] = del[3][3];
                return(dstsq);
            }
        
            // check optimality of face 1-2-4
            d2[9] = del[3][3] - del[3][1];
            d4[9] = del[1][1] - del[3][1];
            e214 = -e124;
            d1[11] = d2[9]*d1[2] + d4[9]*e214;
            d3[14] = d1[11]*d3[4] + d2[11]*e132 + d4[11]*e134;
            if (d1[11]>0.0f && d2[11]>0.0f &&
                d3[14]<=0.0f && d4[11]>0.0f) 
            {
                    *nov = 3;
                    v_index[0][2] = v_index[0][3];
                    v_index[1][2] = v_index[1][3];
                    sum = d1[11] + d2[11] + d4[11];
                    lambda[0] = d1[11]/sum;
                    lambda[1] = d2[11]/sum;
                    lambda[2] = 1.0f - lambda[0] - lambda[1];
                    zsol[0] = y[3][0] + lambda[0]*(y[0][0] - y[3][0]) + lambda[1]*(y[1][0] - y[3][0]);
                    zsol[1] = y[3][1] + lambda[0]*(y[0][1] - y[3][1]) + lambda[1]*(y[1][1] - y[3][1]);
                    zsol[2] = y[3][2] + lambda[0]*(y[0][2] - y[3][2]) + lambda[1]*(y[1][2] - y[3][2]);
                    dstsq = CDYN_DIST_DOT(zsol,zsol);
                    CDYN_DIST_COPY(y[2],y[3]);
                    del[2][0] = del[3][0];
                    del[2][1] = del[3][1];
                    del[2][2] = del[3][3];
                    return(dstsq);
            }

             // check optimality of face 1-3-4
            d3[10] = del[3][3] - del[3][2];
            d4[10] = del[2][2] - del[3][2];
            e314 = -e134;
            d1[12] = d3[10]*d1[4] + d4[10]*e314;
            d2[14] = d1[12]*d2[2] + d3[12]*e123 + d4[12]*e124;
            if (d1[12]>0.0f && d2[14]<=0.0f &&
                d3[12]>0.0f && d4[12]>0.0f) 
            {
                *nov = 3;
                v_index[0][1] = v_index[0][3];
                v_index[1][1] = v_index[1][3];
                sum = d1[12] + d3[12] + d4[12];
                lambda[0] = d1[12]/sum;
                lambda[2] = d3[12]/sum;
                lambda[1] = 1.0f - lambda[0] - lambda[2];
                zsol[0] = y[3][0] + lambda[0]*(y[0][0] - y[3][0]) + lambda[2]*(y[2][0] - y[3][0]);
                zsol[1] = y[3][1] + lambda[0]*(y[0][1] - y[3][1]) + lambda[2]*(y[2][1] - y[3][1]);
                zsol[2] = y[3][2] + lambda[0]*(y[0][2] - y[3][2]) + lambda[2]*(y[2][2] - y[3][2]);
                dstsq = CDYN_DIST_DOT(zsol,zsol);
                CDYN_DIST_COPY(y[1],y[3]);
                del[1][0] = del[3][0];
                del[1][1] = del[3][3];
                del[2][1] = del[3][2];
                return(dstsq);
            }
            
            // check optimality of the hull of all 4 points
            e243 = del[2][1] - del[3][2];
            d4[13] = d2[5]*d4[9] + d3[5]*e243;
            e234 = del[3][1] - del[3][2];
            d3[13] = d2[9]*d3[5] + d4[9]*e234;
            e324 = -e234;
            d2[13] = d3[10]*d2[5] + d4[10]*e324;
            d1[14] = d2[13]*d1[2] + d3[13]*e213 + d4[13]*e214;
            if (d1[14]>0.0f && d2[14]>0.0f &&
                d3[14]>0.0f && d4[14]>0.0f) 
            {
                    sum = d1[14] + d2[14] + d3[14] + d4[14];
                    lambda[0] = d1[14]/sum;
                    lambda[1] = d2[14]/sum;
                    lambda[2] = d3[14]/sum;
                    lambda[3] = 1.0f - lambda[0] - lambda[1] - lambda[2];
                    zsol[0] = lambda[0]*y[0][0] + lambda[1]*y[1][0] + lambda[2]*y[2][0] + lambda[3]*y[3][0];
                    zsol[1] = lambda[0]*y[0][1] + lambda[1]*y[1][1] + lambda[2]*y[2][1] + lambda[3]*y[3][1];
                    zsol[2] = lambda[0]*y[0][2] + lambda[1]*y[1][2] + lambda[2]*y[2][2] + lambda[3]*y[3][2];
                    dstsq = CDYN_DIST_DOT(zsol,zsol);
                    return(dstsq);
            }

            // check optimality of vertex 2
            if (d1[2]<=0.0f && d3[5]<=0.0f && d4[9]<=0.0f) 
            {
                *nov = 1;
                v_index[0][0] = v_index[0][1];
                v_index[1][0] = v_index[1][1];
                lambda[0] = d2[1];
                CDYN_DIST_COPY(zsol,y[1]);
                dstsq = del[1][1];
                CDYN_DIST_COPY(y[0],y[1]);
                del[0][0] = del[1][1];
                return(dstsq);
            }
            
            // check optimality of vertex 3
            if (d1[4]<=0.0f && d2[5]<=0.0f && d4[10]<=0.0f) 
            {
                *nov = 1;
                v_index[0][0] = v_index[0][2];
                v_index[1][0] = v_index[1][2];
                lambda[0] = d3[3];
                CDYN_DIST_COPY(zsol,y[2]);
                dstsq = del[2][2];
                CDYN_DIST_COPY(y[0],y[2]);
                del[0][0] = del[2][2];
                return(dstsq);
            }
        
            // check optimality of vertex 4
            if (d1[8]<=0.0f && d2[9]<=0.0f && d3[10]<=0.0f) 
            {
                *nov = 1;
                v_index[0][0] = v_index[0][3];
                v_index[1][0] = v_index[1][3];
                lambda[0] = d4[7];
                CDYN_DIST_COPY(zsol,y[3]);
                dstsq = del[3][3];
                CDYN_DIST_COPY(y[0],y[3]);
                del[0][0] = del[3][3];
                return(dstsq);
            }
        
            // check optimality of line segment 2-3
            if (d1[6]<=0.0f && d2[5]>0.0f &&
                d3[5]>0.0f && d4[13]<=0.0f) 
            {
                *nov = 2;
                v_index[0][0] = v_index[0][2];
                v_index[1][0] = v_index[1][2];
                sum = d2[5] + d3[5];
                lambda[1] = d2[5]/sum;
                lambda[0] = 1.0f - lambda[1];
                zsol[0] = y[2][0] + lambda[1]*(y[1][0] - y[2][0]);
                zsol[1] = y[2][1] + lambda[1]*(y[1][1] - y[2][1]);
                zsol[2] = y[2][2] + lambda[1]*(y[1][2] - y[2][2]);
                dstsq = CDYN_DIST_DOT(zsol,zsol);
                CDYN_DIST_COPY(y[0],y[2]);
                del[1][0] = del[2][1];
                del[0][0] = del[2][2];
                return(dstsq);
            }
        
            // check optimality of line segment 2-4
            if (d1[11]<=0.0f && d2[9]>0.0f &&
                d3[13]<=0.0f && d4[9]>0.0f) 
            {
                *nov = 2;
                v_index[0][0] = v_index[0][3];
                v_index[1][0] = v_index[1][3];
                sum = d2[9] + d4[9];
                lambda[1] = d2[9]/sum;
                lambda[0] = 1.0f - lambda[1];
                zsol[0] = y[3][0] + lambda[1]*(y[1][0] - y[3][0]);
                zsol[1] = y[3][1] + lambda[1]*(y[1][1] - y[3][1]);
                zsol[2] = y[3][2] + lambda[1]*(y[1][2] - y[3][2]);
                dstsq = CDYN_DIST_DOT(zsol,zsol);
                CDYN_DIST_COPY(y[0],y[3]);
                del[1][0] = del[3][1];
                del[0][0] = del[3][3];
                return(dstsq);
            }
                
            // check optimality of line segment 3-4
            if (d1[12]<=0.0f && d2[13]<=0.0f &&
                d3[10]>0.0f && d4[10]>0.0f) 
            {
                *nov = 2;
                v_index[0][0] = v_index[0][2];
                v_index[0][1] = v_index[0][3];
                v_index[1][0] = v_index[1][2];
                v_index[1][1] = v_index[1][3];
                sum = d3[10] + d4[10];
                lambda[0] = d3[10]/sum;
                lambda[1] = 1.0f - lambda[0];
                zsol[0] = y[3][0] + lambda[0]*(y[2][0] - y[3][0]);
                zsol[1] = y[3][1] + lambda[0]*(y[2][1] - y[3][1]);
                zsol[2] = y[3][2] + lambda[0]*(y[2][2] - y[3][2]);
                dstsq = CDYN_DIST_DOT(zsol,zsol);
                CDYN_DIST_COPY(y[0],y[2]);
                CDYN_DIST_COPY(y[1],y[3]);
                del[0][0] = del[2][2];
                del[1][0] = del[3][2];
                del[1][1] = del[3][3];
                return(dstsq);
            }
            
            // check optimality of face 2-3-4
            if (d1[14]<=0.0f && d2[13]>0.0f &&
                d3[13]>0.0f && d4[13]>0.0f) 
            {
                    *nov = 3;
                    v_index[0][0] = v_index[0][3];
                    v_index[1][0] = v_index[1][3];
                    sum = d2[13] + d3[13] + d4[13];
                    lambda[1] = d2[13]/sum;
                    lambda[2] = d3[13]/sum;
                    lambda[0] = 1.0f - lambda[1] - lambda[2];
                    zsol[0] = y[3][0] + lambda[1]*(y[1][0] - y[3][0]) + lambda[2]*(y[2][0] - y[3][0]);
                    zsol[1] = y[3][1] + lambda[1]*(y[1][1] - y[3][1]) + lambda[2]*(y[2][1] - y[3][1]);
                    zsol[2] = y[3][2] + lambda[1]*(y[1][2] - y[3][2]) + lambda[2]*(y[2][2] - y[3][2]);
                    dstsq = CDYN_DIST_DOT(zsol,zsol);
                    CDYN_DIST_COPY(y[0],y[3]);
                    del[0][0] = del[3][3];
                    del[1][0] = del[3][1];
                    del[2][0] = del[3][2];
                    return(dstsq);
            }
            break;
        } // END case four points

        default: 
        {
            // cDynPrintf(stderr,"Invalid value for nov %d given \n",*nov);
            break;
        }
    } /* END switch */
    return 0; // this line should not be reached.
}


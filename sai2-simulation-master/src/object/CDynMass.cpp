//===========================================================================
/*
    This file is part of the dynamics library.
    Copyright (C) 2014, Artificial Intelligence Laboratory,
    Stanford University. All rights reserved.

    \version   $MAJOR.$MINOR.$RELEASE $Rev: 1242 $
*/
//===========================================================================

//---------------------------------------------------------------------------
#include "utility/CDynLogger.h"
#include "matrix/CDynFrameStack.h"
#include "object/CDynMass.h"
//---------------------------------------------------------------------------

cDynMass::cDynMass()
{
    fs_=NULL;
    p=new Mass;
    p->m_=0.0f;
    p->center_.zero();
    p->inertia_.zero();
    p->rc_=1;
}

//---------------------------------------------------------------------------

cDynMass::~cDynMass()
{
    if (--p->rc_ == 0)
        delete p;
}

//---------------------------------------------------------------------------

cDynMass::cDynMass(cDynMass& mass)
{
    fs_=NULL;
    mass.p->rc_++;
    p=mass.p;
}

//---------------------------------------------------------------------------

cDynMass& cDynMass::operator=(cDynMass& mass)
{
    mass.p->rc_++;
    if (--p->rc_ == 0)
        delete p;
    p=mass.p;

    return *this;
}

//---------------------------------------------------------------------------

// assuming the dimensions stay the same and the body consists of only one
void cDynMass::scale(const double m)
{
    p->inertia_ *= m/p->m_;
    p->m_ = m;
}

//---------------------------------------------------------------------------

void cDynMass::mass(const double m)
{
    cDynVector3 pos;
    if (fs_ == NULL) pos.zero();
    else pos=fs_->top().translation();

    // find inertial term
    double Ixx= m*(pos[1]*pos[1]+pos[2]*pos[2]);
    double Iyy= m*(pos[0]*pos[0]+pos[2]*pos[2]);
    double Izz= m*(pos[0]*pos[0]+pos[1]*pos[1]);
    double Ixy= m*(pos[0]*pos[1]);
    double Ixz= m*(pos[0]*pos[2]);
    double Iyz= m*(pos[1]*pos[2]);
    cDynMatrix3 I;
    I.set(Ixx,	-Ixy,	-Ixz,
         -Ixy,	Iyy,	-Iyz,
         -Ixz,	-Iyz,	Izz);

    p->inertia_ += I;

    // update center of mass
    cDynVector3 sum;
    sum.multiply(p->center_, p->m_);
    cDynVector3 tmpV;
    tmpV.multiply(pos, m);
    sum += tmpV;

    p->m_ += m;
    p->center_.multiply(sum, 1/p->m_);
}

//---------------------------------------------------------------------------

void cDynMass::inertia(const cDynMatrix3& inertia)
{
    if (fs_ == NULL) 
    {
        p->inertia_ += inertia;
    } 
    else 
    {
        cDynMatrix3 R,Rt,tmpM;
        R.set(fs_->top().rotation());
        Rt.transpose(R);
        tmpM.multiply(R, inertia);
        R.multiply(tmpM, Rt);
        // R = R I Rt;
        p->inertia_ += R;
    }
}

//---------------------------------------------------------------------------

void cDynMass::inertia(const cDynVector3& diag)
{
    cDynMatrix3 I;
    I.diagonal(diag);
    inertia(I);
}

//---------------------------------------------------------------------------

void cDynMass::inertia(const double Ixx,const double Iyy,const double Izz)
{
     cDynVector3 I;
     I.set(Ixx,Iyy,Izz);
     inertia(I);
}

//---------------------------------------------------------------------------

void cDynMass::inertia(const double diag[3])
{
    cDynVector3 I;
    I.set(diag[0],diag[1],diag[2]);
    inertia(I);
}

//---------------------------------------------------------------------------

void cDynMass::zero()
{
    // reset all the mass properties of the body making mass, and inertia zero
    p->m_=0.0f;
    p->center_.zero();
    p->inertia_.zero();
}

//---------------------------------------------------------------------------

void cDynMass::cylinder(const double mp,const  double h,const  double r)
{
    //
    // V = (double)M_PI*r*r*h
    //

    double m;
    if (cDynIsDensity(mp)) 
    {
        double p=cDynDensity(mp);
        m=p*(double)(M_PI)*r*r*h;
    } 
    else 
    {
        m=mp;
    }
    
    mass(m);
    double Ixx= (m/12.0f)*(3.0f*r*r + h*h);
    double Izz= 0.5f*m*r*r;
    inertia(Ixx,Ixx,Izz);
}

//---------------------------------------------------------------------------

void cDynMass::cone(const double mp,const  double h,const  double r)
{
    //
    // h in z-axis
    // center of mass h/4 from base in z
    // V = (1.0f/3.0f)*(double)(M_PI)*r*r*h
    //

    double m;
    if (cDynIsDensity(mp)) 
    {
        double p=cDynDensity(mp);
        m=p*(1.0f/3.0f)*(double)(M_PI)*r*r*h;
    } 
    else 
    {
        m=mp;
    }
    
    mass(m);
    double Ixx= (3.0f*m/80.0f)*(4.0f*r*r + h*h);
    double Izz= (3.0f/10.0f)*m*r*r;
    inertia(Ixx,Ixx,Izz);
}

//---------------------------------------------------------------------------


void cDynMass::pyramid(const double mp,const  double a,const  double b,const  double h)
{
    //
    // a in x-axis
    // b in y-axis
    // h in z-axis
    // center of mass h/4 from base in z
    // V = a*b*h/3.0f
    //

    double m;
    if (cDynIsDensity(mp)) 
    {
        double p=cDynDensity(mp);
        m=p*a*b*h/3.0f;
    } 
    else 
    {
        m=mp;
    }
    
    mass(m);
    double Ixx=(m/80.0f)*(4.0f*b*b + 3.0f*h*h);
    double Iyy=(m/80.0f)*(4.0f*a*a + 3.0f*h*h);
    double Izz=(m/20.0f)*(a*a + b*b);
    inertia(Ixx,Iyy,Izz);
}

//---------------------------------------------------------------------------

void cDynMass::block(const double mp,const  double a,const  double b,const  double c)
{
    //
    // a in x-axis
    // b in y-axis
    // c in z-axis
    // V = a*b*c
    //

    double m;
    if (cDynIsDensity(mp)) 
    {
        double p=cDynDensity(mp);
        m=p*a*b*c;
    } 
    else 
    {
        m=mp;
    }
    
    mass(m);
    double m12=m/12.0f;
    double Ixx=m12*(b*b + c*c);
    double Iyy=m12*(a*a + c*c);
    double Izz=m12*(a*a + b*b);
    inertia(Ixx,Iyy,Izz);
}

//---------------------------------------------------------------------------

void cDynMass::sphere(const double mp,const  double r)
{
    //
    // V = (4.0f/3.0f)*(double)(M_PI)*r*r*r
    //

    double m;
    if (cDynIsDensity(mp)) 
    {
        double p=cDynDensity(mp);
        m=p*(4.0f/3.0f)*(double)(M_PI)*r*r*r;
    } 
    else 
    {
        m=mp;
    }
    
    mass(m);
    double I= (2.0f/5.0f)*m*r*r;
    inertia(I,I,I);
}

//---------------------------------------------------------------------------

void cDynMass::hemisphere(const double mp,const  double r)
{
    //
    // center of mass (3/8)*r from base in z
    // V = (2.0f/3.0f)*(double)(M_PI)*r*r*r
    //

    double m;
    if (cDynIsDensity(mp)) 
    {
        double p=cDynDensity(mp);
        m=p*(2.0f/3.0f)*(double)(M_PI)*r*r*r;
    } 
    else 
    {
        m=mp;
    }
    
    mass(m);
    double Ixx=(83.0f/320.0f)*m*r*r;
    double Izz=(2.0f/5.0f)*m*r*r;
    inertia(Ixx,Ixx,Izz);
}

//---------------------------------------------------------------------------


void cDynMass::ellipsoid(const double mp,const  double a,const  double b,const  double c)
{
    //
    // a in x-axis
    // b in y-axis
    // c in z-axis
    // V = (4.0f/3.0f)*(double)(M_PI)*a*b*c
    //
    
    double m;
    if (cDynIsDensity(mp)) 
    {
        double p=cDynDensity(mp);
        m=p*(4.0f/3.0f)*(double)(M_PI)*a*b*c;
    } 
    else 
    {
        m=mp;
    }

    mass(m);
    double m5=m/5.0f;
    double Ixx=m5*(b*b + c*c);
    double Iyy=m5*(a*a + c*c);
    double Izz=m5*(a*a + b*b);
    inertia(Ixx,Iyy,Izz);
}

//---------------------------------------------------------------------------

void cDynMass::rod(const double mp,const  double l)
{
    //
    // l in z-axis
    //

    double m;
    if (cDynIsDensity(mp)) 
    {
        double p=cDynDensity(mp);
        m=p*l;
    } 
    else 
    {
        m=mp;
    }
    
    mass(m);
    double I=m*l*l/12.0f;
    inertia(I,I,0.0f);
}

//---------------------------------------------------------------------------

void cDynMass::disk(const double mp,const double r)
{
    //
    // A = M_PI*r*r
    //

    double m;
    if (cDynIsDensity(mp)) 
    {
        double p=cDynDensity(mp);
        m=p*(double)(M_PI)*r*r;
    } 
    else 
    {
        m=mp;
    }
    
    mass(m);
    double Ixx=0.25f*m*r*r;
    double Izz=0.5f*m*r*r;
    inertia(Ixx,Ixx,Izz);
}

//---------------------------------------------------------------------------

void cDynMass::plate(const double mp,const  double a,const  double b)
{
    //
    // a in x-axis
    // b in y-axis
    // A = a*b
    //

    double m;
    if (cDynIsDensity(mp)) 
    {
        double p=cDynDensity(mp);
        m=p*a*b;
    } 
    else 
    {
        m=mp;
    }
    
    mass(m);
    double m12=m/12.0f;
    double Ixx=m12*b*b;
    double Iyy=m12*a*a;
    double Izz=m12*(a*a + b*b);
    inertia(Ixx,Iyy,Izz);
}

//---------------------------------------------------------------------------

void cDynMass::cylinderShell(const double mp,const  double h,const  double r)
{
    //
    // h in z-axis
    // A = 2.0f*(double)(M_PI)*r*h
    //

    double m;
    if (cDynIsDensity(mp)) 
    {
        double p=cDynDensity(mp);
        m=p*2.0f*(double)(M_PI)*r*h;
    } 
    else 
    {
        m=mp;
    }
    
    mass(m);
    double Ixx=(m/12.0f)*(6*r*r + h*h);
    double Izz=m*r*r;
    inertia(Ixx,Ixx,Izz);
}

//---------------------------------------------------------------------------

void cDynMass::coneShell(const double mp,const  double h,const  double r)
{
    //
    // h in z-axis
    // center of mass (1/3)*h from base in z
    // A = (double)(M_PI)*r*cDynSqrt(r*r + h*h)
    //
    
    double m;
    if (cDynIsDensity(mp)) 
    {
        double p=cDynDensity(mp);
        m=p*(double)(M_PI)*r*cDynSqrt(r*r + h*h);
    } 
    else 
    {
        m=mp;
    }
    
    mass(m);
    double Ixx=(m/18.0f)*(4.5f*r*r + h*h);
    double Izz=0.5f*m*r*r;
    inertia(Ixx,Ixx,Izz);
}

//---------------------------------------------------------------------------

void cDynMass::sphereShell(const double mp,const  double r)
{
    //
    // A = 4.0f*(double)(M_PI)*r*r
    //

    double m;
    if (cDynIsDensity(mp)) 
    {
        double p=cDynDensity(mp);
        m=p*4.0f*(double)(M_PI)*r*r;
    } 
    else 
    {
        m=mp;
    }
    
    mass(m);
    double I=(2.0f/3.0f)*m*r*r;
    inertia(I,I,I);
}

//---------------------------------------------------------------------------

void cDynMass::hemisphereShell(const double mp,const  double r)
{
    //
    // center of mass r/2.0f from base in z
    // A = 2.0f*(double)(M_PI)*r*r
    //

    double m;
    if (cDynIsDensity(mp)) 
    {
        double p=cDynDensity(mp);
        m=p*2.0f*(double)(M_PI)*r*r;
    } 
    else 
    {
        m=mp;
    }
    
    mass(m);
    double Ixx=(5.0f/12.0f)*m*r*r;
    double Izz=(2.0f/3.0f)*m*r*r;
    inertia(Ixx,Ixx,Izz);
}

//---------------------------------------------------------------------------
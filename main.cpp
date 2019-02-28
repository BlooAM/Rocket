#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <windows.h>
#include <vector>
#include "winbgi2.h"
#ifndef M_PI
const double M_PI=3.14159265358979323846;
#endif

using namespace std;

//Define global rocket properties
const double hTriangle = 10, hRectangle = 20, hTrapezoid = 10, b = 20, a = 10;
double xc = 0, yc = ((a+b)*(hTrapezoid/2)*hTrapezoid*(b+2*a)/(3*(b+a)) + a*hRectangle*(hTrapezoid+hRectangle/2)+a*(hTriangle/2)*(hTrapezoid+hRectangle+hTriangle/3))/((a+b)*hTrapezoid/2 + a*hRectangle + a*hTriangle/2);
double mass = 1, J = 10, Fg = 2, c=50 ;   //J = 1000
double xDest=250+c,yDest=250-c;
double dt = 0.01;
int GraphicsON=0;

//Define graphics window properties
const double width = 500, height = 500;


//Declaration of all functions
void integrateTimeStep(double *,double *,double*,double*);
void integrateTimeStep_b(double *, double *, double *, double *, double *,double *, double *, double *);
void ObjectiveFunction(double *,double *,int);
void ObjectiveFunction_b(double *, double *, double *, double *,int);
void drawRocket(double *,double *);

int main()
{
    int nTimeSteps = 50,max_iter=100000, max_inner_iter = 10; //max_iter = 100;
    double initState[] = {width/2,height/2,0,0,M_PI/2,0}, *OF = new double [nTimeSteps], *OFb = new double [nTimeSteps];
    double Objb=1,Obj = 0, alfa = 1;
    double **v = new double *[nTimeSteps], **vb = new double *[nTimeSteps], **force = new double *[nTimeSteps-1], **forceb = new double *[nTimeSteps-1], **v_s = new double *[2];
    double temp_f[2],temp_obj;
    GraphicsON = 1;


    for(int i=0; i<nTimeSteps; i++)
    {
        v[i] = new double [6];
        vb[i] = new double [6];
        force[i] = new double [2];
        forceb[i] = new double [2];
    }
    v_s[0] = new double [6];
    v_s[1] = new double [6];
    for(int i=0; i<2; i++)
    {
        v_s[i][0] = 250;
        v_s[i][1] = 250;
        v_s[i][2] = 0;
        v_s[i][3] = 0;
        v_s[i][4] = M_PI/2;
        v_s[i][5] = 0;
    }
    for(int i=0; i<nTimeSteps; i++)
    {
        OF[i] = 0;
        OFb[i] = 0;
        for(int k=0; k<6; k++)
        {
            v[i][k] = 0;
            vb[i][k] = 0;
            if(k<2)
            {
                forceb[i][k] = 0;
                force[i][k] = 1;
            }
        }
    }
    v[0][0] = v_s[0][0];
    v[0][1] = v_s[0][1];
    v[0][2] = v_s[0][2];
    v[0][3] = v_s[0][3];
    v[0][4] = v_s[0][4];
    v[0][5] = v_s[0][5];
    if(GraphicsON==1 || GraphicsON==2)
    {
        graphics(width, height);
        drawRocket(v[0],force[0]);
        circle(xDest,height-yDest,3);
        //Sleep(100);
        clear();
    }

    for(int iter = 0; iter<max_iter; iter++)
    {
        temp_obj = 1e100;
        for(int inner_iter = 0; inner_iter<max_inner_iter; inner_iter++)
        {
            //FORWARD
            Obj = 0;
            OF[0] = ((v[0][0]-xDest)*(v[0][0]-xDest)+(v[0][1]-yDest)*(v[0][1]-yDest));
            for(int i=1; i<nTimeSteps; i++)
            {
                integrateTimeStep(v[i-1],v[i],force[i-1],&OF[i]);
                if(GraphicsON==0)   cout<<v[i][0]<<"\t"<<v[i][1]<<"\t"<<v[i][2]<<"\t"<<v[i][3]<<endl;
            }

            ObjectiveFunction(OF, &Obj, nTimeSteps);
            if(Obj<temp_obj)
            {
                temp_f[0] = force[1][0];
                temp_f[1] = force[1][1];
            }
            for(int i=0; i<nTimeSteps; i++)
            {
                OFb[i] = 0;
                for(int k=0; k<6; k++)
                {
                    vb[i][k] = 0;
                    if(k<2)
                    {
                        forceb[i][k] = 0;
                    }
                }
            }
            if(GraphicsON==2)
            {
                drawRocket(v[0],force[0]);
                circle(xDest,height-yDest,3);
                for(int i=1; i<nTimeSteps; i++)
                {
                    line(v[i-1][0],height-v[i-1][1],v[i][0],height-v[i][1]);
                    //cout<<"\t"<<v[i-1][0]<<"\t"<<v[i-1][1]<<"\t"<<v[i][0]<<"\t"<<v[i][1]<<endl;
                    //cout<<"Iter: "<<inner_iter<<"\tForces: "<<force[i][0]<<"\t"<<force[i][1]<<endl;
                    //Sleep(100);
                    //clear();
                }
                cout<<"Obj funct: "<<Obj<<endl<<endl;
                //wait();
            }

            //ADJOINT
            ObjectiveFunction_b(OF,OFb, &Obj, &Objb, nTimeSteps);
            for(int i=nTimeSteps-2; i>-1; --i)
            {
                integrateTimeStep_b(v[i],vb[i],v[i+1],vb[i+1],force[i],forceb[i+1],&OF[i],&OFb[i]);
            }
            for(int i=0; i<nTimeSteps; i++)
            {
                force[i][0] = force[i][0] - alfa*forceb[i][0]; //force[i][0] = force[i][0] - alfa*forceb[i][0];
                force[i][1] = force[i][1] - alfa*forceb[i][1]; //force[i][1] = force[i][1] - alfa*forceb[i][1];
            }
            //if(GraphicsON==0)   cout<<"Init state: "<<v[0][0]<<"\t"<<v[0][1]<<endl;
            if(GraphicsON==0)  system("pause");
        }
        if(GraphicsON==0)   cout<<"AND THE WINNER IS: "<<temp_f[0]<<"\t"<<temp_f[1]<<endl;
        if(GraphicsON==2)
        {
            Sleep(1000);
            wait();
            clear();
        }
        for(int i=0; i<1; i++)
        {
            v_s[0][0] = v_s[1][0];
            v_s[0][1] = v_s[1][1];
            v_s[0][2] = v_s[1][2];
            v_s[0][3] = v_s[1][3];
            v_s[0][4] = v_s[1][4];
            v_s[0][5] = v_s[1][5];
            integrateTimeStep(v_s[0],v_s[1],temp_f,&OF[0]);
            if(GraphicsON==1)
            {
                drawRocket(v_s[1],temp_f);
                circle(xDest,height-yDest,3);
                //Sleep(100);
                clear();
            }
        }
        v[0][0] = v_s[1][0];
        v[0][1] = v_s[1][1];
        v[0][2] = v_s[1][2];
        v[0][3] = v_s[1][3];
        v[0][4] = v_s[1][4];
        v[0][5] = v_s[1][5];

        if(GraphicsON==0)   cout<<"Forces: "<<temp_f[0]<<"\t"<<temp_f[1]<<"\t Diff: "<<temp_f[0]-temp_f[1]<<"\tSum: "<<temp_f[0]+temp_f[1]-Fg<<endl;
        if(GraphicsON==0)   cout<<"x = "<<v[0][0]<<"\t y = "<<v[0][1]<<endl;
        //if(iter%1==0) system("pause");
    }

    for (int i=0; i<nTimeSteps-1; i++)    delete [] v[i],vb[i],force[i],forceb[i];
    for (int i=0; i<2; i++) delete [] v_s[i];
    delete [] v[nTimeSteps-1], vb[nTimeSteps-1];
    delete [] v,vb,force,forceb,v_s;
    delete [] OF, OFb;
    if(GraphicsON==0)   system("pause");
    return 0;
}
//******************************************************************************
void integrateTimeStep(double *v1,double *v2,double *force,double *OF)
{
    double F = force[0]+force[1], torque = (force[1]-force[0])*(a/2);
    double Fx, Fy;
    Fx = F*cos(v1[4]);
    Fy = F*sin(v1[4]) - Fg;
    v2[2] = v1[2]+(Fx/mass)*dt;
    v2[3] = v1[3]+(Fy/mass)*dt;
    v2[0] = v1[0]+v2[2]*dt; //Semi-implicit
    v2[1] = v1[1]+v2[3]*dt; //Euler method !
    v2[5] = v1[5]+(torque/J)*dt;
    v2[4] = v1[4]+v2[5]*dt;
    *OF = ((v2[0]-xDest)*(v2[0]-xDest)+(v2[1]-yDest)*(v2[1]-yDest));
}
//******************************************************************************
void integrateTimeStep_b(double *v1, double *v1b, double *v2, double *v2b, double *force, double *forceb, double *OF, double *OFb)
{
    double F = force[0] + force[1];
    double Fb = 0.0;
    double torque = (force[1]-force[0])*(a/2);
    double torqueb = 0.0;
    double Fx, Fy;
    double Fxb, Fyb;
    double tempb;
    Fx = F*cos(v1[4]);
    Fy = F*sin(v1[4]) - Fg;
    v2[2] = v1[2] + Fx/mass*dt;
    v2[3] = v1[3] + Fy/mass*dt;
    v2[0] = v1[0] + v2[2]*dt;
    v2[1] = v1[1] + v2[3]*dt;
    v2[5] = v1[5] + torque/J*dt;
    v2[4] = v1[4] + v2[5]*dt;
    v2b[0] = v2b[0] + 2*(v2[0]-xDest)*(*OFb);
    v2b[1] = v2b[1] + 2*(v2[1]-yDest)*(*OFb);
    *OFb = 0.0;
    v1b[4] = v1b[4] + v2b[4];
    v2b[5] = v2b[5] + dt*v2b[4];
    v2b[4] = 0.0;
    v1b[5] = v1b[5] + v2b[5];
    torqueb = dt*v2b[5]/J;
    v2b[5] = 0.0;
    v1b[1] = v1b[1] + v2b[1];
    v2b[3] = v2b[3] + dt*v2b[1];
    v2b[1] = 0.0;
    v1b[0] = v1b[0] + v2b[0];
    v2b[2] = v2b[2] + dt*v2b[0];
    v2b[0] = 0.0;
    v1b[3] = v1b[3] + v2b[3];
    Fyb = dt*v2b[3]/mass;
    v2b[3] = 0.0;
    v1b[2] = v1b[2] + v2b[2];
    Fxb = dt*v2b[2]/mass;
    v2b[2] = 0.0;
    Fb = cos(v1[4])*Fxb + sin(v1[4])*Fyb;
    v1b[4] = v1b[4] + F*cos(v1[4])*Fyb;
    v1b[4] = v1b[4] - F*sin(v1[4])*Fxb;
    tempb = a*torqueb/2;
    forceb[1] = forceb[1] + tempb;
    forceb[0] = forceb[0] - tempb;
    forceb[0] = forceb[0] + Fb;
    forceb[1] = forceb[1] + Fb;
}
//******************************************************************************
void ObjectiveFunction(double *OF,double *OBJ,int N)
{
    for(int i=0; i<N; i++)  *OBJ = *OBJ + OF[i];
}
//******************************************************************************
void ObjectiveFunction_b(double *OF, double *OFb, double *OBJ, double *OBJb,int N)
{
    for (int i = N-1; i > -1; --i)
        OFb[i] = OFb[i] + *OBJb;
}
//******************************************************************************
void drawRocket(double *v,double *force)
{
    double H = hTrapezoid+hRectangle+hTriangle, fi1 = atan((a/2)/hTriangle), l1 = sqrt((a/2)*(a/2)+hTriangle*hTriangle);
    double l2 = sqrt(((b-a)*(b-a)/4)+hTrapezoid*hTrapezoid), fi2 = atan(hTrapezoid/((b-a)/2));
    double r = 2, amp = 5;
    double x = v[0], y = height-v[1];
    double xTop = x+(H-yc)*cos(v[4]), yTop = y-(H-yc)*sin(v[4]);
    double xTopR = xTop+l1*cos(M_PI-fi1-v[4]), yTopR = yTop+l1*sin(M_PI-fi1-v[4]);
    double xTopL = xTop-l1*cos(v[4]-fi1), yTopL = yTop+l1*sin(v[4]-fi1);
    double xMiddleR = xTopR-hRectangle*cos(v[4]), yMiddleR = yTopR+hRectangle*sin(v[4]);
    double xMiddleL = xTopL-hRectangle*cos(v[4]), yMiddleL = yTopL+hRectangle*sin(v[4]);
    double xBottomR = xMiddleR+l2*cos(fi2-v[4]+M_PI/2), yBottomR = yMiddleR + l2*sin(fi2-v[4]+M_PI/2);
    double xBottomL = xMiddleL-l2*cos(v[4]+fi2-M_PI/2), yBottomL = yMiddleL + l2*sin(v[4]+fi2-M_PI/2);
    /*cout<<"xc = "<<x<<"\tyc = "<<y<<endl;
    cout<<"xtop = "<<xTop<<"\tytop = "<<yTop<<endl;
    cout<<"xtopr = "<<xTopR<<"\tytopr = "<<yTopR<<endl;
    cout<<"xtopl = "<<xTopL<<"\tytopl = "<<yTopL<<endl;
    cout<<"xmiddler = "<<xMiddleR<<"\tymiddler = "<<yMiddleR<<endl;
    cout<<"xmiddlel = "<<xMiddleL<<"\tymiddlel = "<<yMiddleL<<endl;
    cout<<"xbottomr = "<<xBottomR<<"\tybottomr = "<<yBottomR<<endl;
    cout<<"xbottoml = "<<xBottomL<<"\tybottoml = "<<yBottomL<<endl;*/

    circle(x,y,r);
    /*circle(xTop,yTop,r);
    circle(xTopR,yTopR,r);
    circle(xTopL,yTopL,r);
    circle(xMiddleR,yMiddleR,r);
    circle(xMiddleL,yMiddleL,r);
    circle(xBottomR,yBottomR,r);
    circle(xBottomL,yBottomL,r);*/
    line(xTop,yTop,xTopR,yTopR);
    line(xTopR,yTopR,xMiddleR,yMiddleR);
    line(xMiddleR,yMiddleR,xBottomR,yBottomR);
    line(xBottomR,yBottomR,xBottomL,yBottomL);
    line(xBottomL,yBottomL,xMiddleL,yMiddleL);
    line(xMiddleL,yMiddleL,xTopL,yTopL);
    line(xTopL,yTopL,xTop,yTop);
    line(xBottomL+0.25*(xBottomR-xBottomL)-amp*force[0]*cos(v[4]),yBottomL+0.25*(yBottomR-yBottomL)+amp*force[0]*sin(v[4]),xBottomL+0.25*(xBottomR-xBottomL),yBottomL+0.25*(yBottomR-yBottomL));
    line(xBottomL+0.75*(xBottomR-xBottomL)-amp*force[1]*cos(v[4]),yBottomL+0.75*(yBottomR-yBottomL)+amp*force[1]*sin(v[4]),xBottomL+0.75*(xBottomR-xBottomL),yBottomL+0.75*(yBottomR-yBottomL));

}

void integrateTimeStep(double *v1,double *v2,double *force,double *OF)
{
	double hTriangle = 10, hRectangle = 20, hTrapezoid = 10, b = 20, a = 10;
	double xc = 0, yc = ((a+b)*(hTrapezoid/2)*hTrapezoid*(b+2*a)/(3*(b+a)) + a*hRectangle*(hTrapezoid+hRectangle/2)+a*(hTriangle/2)*(hTrapezoid+hRectangle+hTriangle/3))/((a+b)*hTrapezoid/2 + a*hRectangle + a*hTriangle/2);
	double mass = 1 , J = 1, Fg = 2 ;  //J = 1000
	double xDest=500/2+100,yDest=500/2+100;

	double width = 500, height = 500;
    double dt = 0.01;
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

void integrateTimeStep_b(double *v1, double *v2, double *v2b, double *force, double *forceb, double *OF, double *OFb) 
{
    double dt = 0.01;
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
    v2b[5] = v2b[5] + dt*v2b[4];
    v2b[4] = 0.0;
    torqueb = dt*v2b[5]/J;
    v2b[5] = 0.0;
    v2b[3] = v2b[3] + dt*v2b[1];
    v2b[1] = 0.0;
    v2b[2] = v2b[2] + dt*v2b[0];
    v2b[0] = 0.0;
    Fyb = dt*v2b[3]/mass;
    v2b[3] = 0.0;
    Fxb = dt*v2b[2]/mass;
    Fb = cos(v1[4])*Fxb + sin(v1[4])*Fyb;
    tempb = a*torqueb/2;
    forceb[1] = forceb[1] + tempb;
    forceb[0] = forceb[0] - tempb;
    forceb[0] = forceb[0] + Fb;
    forceb[1] = forceb[1] + Fb;
}
void integrateTimeSteps(double **v,double *force,double *OF, double *OBJ,int N)
{
	double hTriangle = 10, hRectangle = 20, hTrapezoid = 10, b = 20, a = 10;
	double xc = 0, yc = ((a+b)*(hTrapezoid/2)*hTrapezoid*(b+2*a)/(3*(b+a)) + a*hRectangle*(hTrapezoid+hRectangle/2)+a*(hTriangle/2)*(hTrapezoid+hRectangle+hTriangle/3))/((a+b)*hTrapezoid/2 + a*hRectangle + a*hTriangle/2);
	double mass = 1 , J = 1, Fg = 2 ;  
	double xDest=500/2+100,yDest=500/2+100;
	double width = 500, height = 500;

    double dt = 0.1;
    double F = force[0]+force[1], torque = (force[1]-force[0])*(a/2);
    double Fx, Fy;
    OF[0] =((v[0][0]-xDest)*(v[0][0]-xDest)+((height-v[0][1])-yDest)*((height-v[0][1])-yDest));
    *OBJ = OF[0];
    for(int i=1; i<N; i++)
    {
        Fx = F*cos(v[i-1][4]);
        Fy = F*sin(v[i-1][4]) - Fg;
        v[i][2] = v[i-1][2]+(Fx/mass)*dt;
        v[i][3] = v[i-1][3]+(Fy/mass)*dt;
        v[i][0] = v[i-1][0]+v[i][2]*dt; 
        v[i][1] = v[i-1][1]+v[i][3]*dt; 
        v[i][5] = v[i-1][5]+(torque/J)*dt;
        v[i][4] = v[i-1][4]+v[i][5]*dt;
        OF[i] = ((v[i][0]-xDest)*(v[i][0]-xDest)+((height-v[i][1])-yDest)*((height-v[i][1])-yDest));
        *OBJ = *OBJ + OF[i];
    }
}


void integrateTimeSteps_b(double **v, double **vb, double *force, double *forceb, double *OF, double *OFb, double *OBJ, double *OBJb, int N) 
{
    double dt = 0.1;
    double F = force[0] + force[1];
    double Fb = 0.0;
    double torque = (force[1]-force[0])*(a/2);
    double torqueb = 0.0;
    double Fx, Fy;
    double Fxb, Fyb;
    double tempb;
    {
      double tmp;
      double tmp0;
      double tmp1;
      double tmp2;
      double tmp3;
      double tmp4;
      double tmpb;
      double tmpb0;
      double tmpb1;
      double tmpb2;
      double tmpb3;
      double tmpb4;
      for (int i = 1; i < N; ++i) {
          Fx = F*cos(v[i-1][4]);
          Fy = F*sin(v[i-1][4]) - Fg;
          tmp = v[i - 1][2] + Fx/mass*dt;
          pushreal8(v[i][2]);
          v[i][2] = tmp;
          tmp0 = v[i - 1][3] + Fy/mass*dt;
          pushreal8(v[i][3]);
          v[i][3] = tmp0;
          tmp1 = v[i - 1][0] + v[i][2]*dt;
          pushreal8(v[i][0]);
          v[i][0] = tmp1;
          tmp2 = v[i - 1][1] + v[i][3]*dt;
          pushreal8(v[i][1]);
          v[i][1] = tmp2;
          tmp3 = v[i - 1][5] + torque/J*dt;
          pushreal8(v[i][5]);
          v[i][5] = tmp3;
          tmp4 = v[i - 1][4] + v[i][5]*dt;
          pushreal8(v[i][4]);
          v[i][4] = tmp4;
      }
      **vb = 0.0;
      *OFb = 0.0;
      torqueb = 0.0;
      Fb = 0.0;
      for (int i = N-1; i > 0; --i) {
          OFb[i] = OFb[i] + *OBJb;
          vb[i][0] = vb[i][0] + 2*(v[i][0]-xDest)*OFb[i];
          vb[i][1] = vb[i][1] - 2*(height-yDest-v[i][1])*OFb[i];
          OFb[i] = 0.0;
          popreal8(&v[i][4]);
          tmpb = vb[i][4];
          vb[i][4] = 0.0;
          vb[i - 1][4] = vb[i - 1][4] + tmpb;
          vb[i][5] = vb[i][5] + dt*tmpb;
          popreal8(&v[i][5]);
          tmpb0 = vb[i][5];
          vb[i][5] = 0.0;
          vb[i - 1][5] = vb[i - 1][5] + tmpb0;
          torqueb = torqueb + dt*tmpb0/J;
          popreal8(&v[i][1]);
          tmpb1 = vb[i][1];
          vb[i][1] = 0.0;
          vb[i - 1][1] = vb[i - 1][1] + tmpb1;
          vb[i][3] = vb[i][3] + dt*tmpb1;
          popreal8(&v[i][0]);
          tmpb2 = vb[i][0];
          vb[i][0] = 0.0;
          vb[i - 1][0] = vb[i - 1][0] + tmpb2;
          vb[i][2] = vb[i][2] + dt*tmpb2;
          popreal8(&v[i][3]);
          tmpb3 = vb[i][3];
          vb[i][3] = 0.0;
          vb[i - 1][3] = vb[i - 1][3] + tmpb3;
          Fyb = dt*tmpb3/mass;
          popreal8(&v[i][2]);
          tmpb4 = vb[i][2];
          vb[i][2] = 0.0;
          vb[i - 1][2] = vb[i - 1][2] + tmpb4;
          Fxb = dt*tmpb4/mass;
          Fb = Fb + cos(v[i-1][4])*Fxb + sin(v[i-1][4])*Fyb;
          vb[i - 1][4] = vb[i - 1][4] + F*cos(v[i-1][4])*Fyb;
          vb[i - 1][4] = vb[i - 1][4] - F*sin(v[i-1][4])*Fxb;
      }
    }
    *OBJb = 0.0;
    tempb = a*torqueb/2;
    forceb[1] = forceb[1] + tempb;
    forceb[0] = forceb[0] - tempb;
    forceb[0] = forceb[0] + Fb;
    forceb[1] = forceb[1] + Fb;
}
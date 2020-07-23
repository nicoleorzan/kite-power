#include <math.h>

#ifndef WIND
#define WIND

#define PI 3.1415926535897932384626433

#define epsilon_wind 0.1
#define epsilon_wind_hard 0.2
#define Lx 50.0
#define Ly 50.0
#define k_wind 0.5

#define mean_wind 25
#define tau 5.0
#define sigma 1.0

//from Houska - Nonlinear MPC of kites under varying wind conditions for a new class of large scale wind power generators
//#define w0 1//10.0        // m/s
#define h0 100              // m
#define hr 0.1              // m, roughness length

#define hh 10

void streamfunction2d(double *rk, double *W){ // kp kiteposition in (x,z)

    W[0] = 0.5*k_wind*rk[1]*(2*epsilon_wind*sin(PI*rk[0]/Lx)*sin(PI*rk[1]/Ly) + \
        epsilon_wind*PI*rk[1]/Ly*sin(PI*rk[0]/Lx)*cos(PI*rk[1]/Ly) + 2);
    W[1] = -k_wind*epsilon_wind*PI*rk[1]*rk[1]/(2*Lx)*cos(PI*rk[0]/Lx)*sin(PI*rk[1]/Ly);

}

void streamfunction2d_hard(double *rk, double *W){ // kp kiteposition in (x,z)

    W[0] = 0.5*k_wind*rk[1]*(2*epsilon_wind_hard*sin(PI*rk[0]/Lx)*sin(PI*rk[1]/Ly) + \
        epsilon_wind_hard*PI*rk[1]/Ly*sin(PI*rk[0]/Lx)*cos(PI*rk[1]/Ly) + 2);
    W[1] = -k_wind*epsilon_wind_hard*PI*rk[1]*rk[1]/(2*Lx)*cos(PI*rk[0]/Lx)*sin(PI*rk[1]/Ly);

}

void streamfunction3d_hard(double *rk, double *W){ // kp kiteposition in (x,y,z)

    W[0] = 0.5*k_wind*rk[2]*(2*epsilon_wind_hard*sin(PI*rk[0]/Lx)*sin(PI*rk[2]/Ly) + \
      epsilon_wind_hard*PI*rk[2]/Ly*sin(PI*rk[0]/Lx)*cos(PI*rk[2]/Ly) + 2);
    W[1] = 0;
    W[2] = -k_wind*epsilon_wind_hard*PI*rk[2]*rk[2]/(2*Lx)*cos(PI*rk[0]/Lx)*sin(PI*rk[2]/Ly);
}

double ornstein_uhlenbeck(double * fluct){

    double u1 = (double)rand() / (double)RAND_MAX;
    double u2 = (double)rand() / (double)RAND_MAX;
    double z1 = sqrt(-2.0*log(u1))*cos(2*PI*u2);
    //double z2 = sqrt(-2.0*log(u1))*sin(2*PI*u2);

    double der_fluct = -*fluct/tau*hh + sigma*z1*sqrt(hh); // ????
    *fluct = *fluct + der_fluct;

    return mean_wind + *fluct;
}

double wind_shear(double z){
    return mean_wind*log(z/hr)/log(h0/hr);
}

#endif

#include <math.h>

#ifndef WIND
#define WIND

#define PI 3.1415926535897932384626433

#define epsilon_wind 0.1
#define Lx 50
#define Ly 50
#define k 0.5

#define mean_wind 25
#define tau 5.0
#define sigma 1.0

//from Houska - Nonlinear MPC of kites under varying wind conditions for a new class of large scale wind power generators
//#define w0 1//10.0        // m/s
#define h0 100              // m
#define hr 0.1              // m, roughness length

void streamfunction(double *rk, double *ux, double *uy){ // kp kiteposition in (x,z)

    *ux = k*rk[1]*(1+epsilon_wind*sin(PI*rk[0]/Lx)*cos(PI*rk[1]/Ly)*PI/Ly);
    *uy = -(k*rk[1]*rk[1]/2*epsilon_wind*sin(PI*rk[1]/Ly)*cos(PI*rk[0]/Lx));
}

double ornstein_uhlenbeck(double * fluct){

    double u1 = (double)rand() / (double)RAND_MAX;
    double u2 = (double)rand() / (double)RAND_MAX;
    double z1 = sqrt(-2.0*log(u1))*cos(2*PI*u2);
    //double z2 = sqrt(-2.0*log(u1))*sin(2*PI*u2);

    double der_fluct = -*fluct/tau*h + sigma*z1*sqrt(h); // ????
    *fluct = *fluct + der_fluct;

    return mean_wind + *fluct;
}

double wind_shear(double z){
    return mean_wind*log(z/hr)/log(h0/hr);
}

#endif

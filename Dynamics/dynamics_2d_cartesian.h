#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "constants.h"

#ifndef __dynamics__
#define __dynamics__

#define PI 3.1415926535897932384626433

//#define DEBUG

#define h 0.0001                // Integration step
#define coeff_friction 0.4       // Block friction coefficient
#define m_block 1000.0          // Block mass

double W;
double Va_mod;
double T_denom;
double denom1, denom2;
double T1, T2;
double va[2];
double L[2], Lc;
double D[2], Dc;
double t2_mod, t3_mod;
double Fg[2] = {0, -m*g};
double F_aer[2];
double F_attrito[2];
double F_attrito_mod;
double N;
double beta;
double Tension[2];
double r_diff_mod;

double T_orizz;
double F_attr_sign;
double Ftot[3]; 

void variables_initialization(double * rk, double * vk, double * ak, double theta,
                              double * r_block, double * v_block, double * a_block){
    r_block[0] = 0;
    r_block[1] = 0;

    v_block[0] = 0;
    v_block[1] = 0;

    a_block[0] = 0;
    a_block[1] = 0;

    rk[0] = r_block[0] + R*sin(theta);
    rk[1] = r_block[1] + R*cos(theta);

    vk[0] = 0;
    vk[1] = 0;

    ak[0] = 0;
    ak[1] = 0;
}

void integration_trajectory(double * rk, double * vk, double * ak, // Kite variables
                            double * r_block, double * v_block, double * a_block, // Block variables
                            double * r_diff, double * v_diff, double * a_diff, 
                            double * theta,
                            int alpha,
                            double * W, double * lift, double * drag,
                            double * T, int it){

r_diff[0] = rk[0] - r_block[0];
r_diff[1] = rk[1] - r_block[1];

*theta = atan(r_diff[1]/r_diff[0]);

r_diff_mod = sqrt(r_diff[0]*r_diff[0] + r_diff[1]*r_diff[1]);
                       
va[0] = vk[0] - W[0];              // Apparent velocity on x
va[1] = vk[1] - W[1];              // Apparent velocity on z

Va_mod = sqrt(va[0]*va[0] + va[1]*va[1]);

// Computing Lift and Drag     

beta = atan2(va[1], va[0]);

Lc = 0.5*rho*CL_alpha[alpha]*A*Va_mod*Va_mod;
*lift = Lc;

Dc = 0.5*rho*CD_alpha[alpha]*A*Va_mod*Va_mod;
*drag = Dc;

// Lift (x, z)

if (beta > - PI/2. && beta < PI/2.){
    L[0] = Lc*cos(beta + PI/2.);
    L[1] = Lc*sin(beta + PI/2.);
} else {
    L[0] = Lc*cos(beta - PI/2.);
    L[1] = Lc*sin(beta - PI/2.);
}

// Drag (x, z)

D[0] = Dc*cos(beta + PI);
D[1] = Dc*sin(beta + PI);       

F_aer[0] = L[0] + D[0];
F_aer[1] = L[1] + D[1];

if (v_block[0] != 0){

    denom1 = R*(m+m_block)/(m*m_block)
            - cos(*theta)/m_block*(rk[1] - coeff_friction*r_diff[0]*v_block[0]/fabs(v_block[0]));

    T1 = (F_aer[0]*r_diff[0] + F_aer[1]*r_diff[1])/m
        + (v_diff[0]*v_diff[0] + v_diff[1]*v_diff[1]) 
        - g*(r_diff[1] + coeff_friction*r_diff[0]*v_block[0]/fabs(v_block[0]));

    T1 = T1/denom1;

    *T = T1;

    Tension[0] = *T*sin(*theta);
    Tension[1] = *T*cos(*theta);
 
    N = m_block*g - Tension[1];

    F_attrito_mod = coeff_friction*fabs(N);

    F_attrito[0] = -F_attrito_mod*v_block[0]/fabs(v_block[0]*v_block[0] + v_block[1]*v_block[1]);

    //printf("vblock0=%f\n", v_block[0]); 
    //printf("denom1=%f\n", denom1); 
}
else if (v_block[0] == 0){
    
    denom2 = R*(m+m_block)/(m*m_block)
        - sin(*theta)/m_block*r_diff[0] - cos(*theta)*rk[1]/m_block;
    T2 = (F_aer[0]*r_diff[0] + F_aer[1]*r_diff[1])/m
        + (v_diff[0]*v_diff[0] + v_diff[1]*v_diff[1]) 
        - g*r_diff[1];

    T2 = T2/denom2;

    *T = T2;

    Tension[0] = *T*sin(*theta);
    Tension[1] = *T*cos(*theta);
 
    N = m_block*g - Tension[1];

    F_attrito[0] = -Tension[0];
    //printf("vblock0=%f\n", v_block[0]); 
}

// sulla forza d'attrito va messo il meno a mano?? forse si!

a_block[0] = (Tension[0] + F_attrito[0])/m_block;
a_block[1] = 0;

v_block[0] = v_block[0] + h*a_block[0]; 
v_block[1] = v_block[1] + h*a_block[1];

if (it%10000 == 0){
    //printf("Lx=%f, Ly=%f, Lz=%f\n", L[0], L[1], L[2]);
    //printf("Tx=%f, Ty=%f, Tz=%f\n", Tension[0], Tension[1], Tension[2]);
    //printf("Fx=%f, vblockx=%f\n", F_attrito[0], v_block[0]);
}

r_block[0] = r_block[0] + h*v_block[0]; 
r_block[1] = r_block[1] + h*v_block[1];

ak[0] = (D[0] + L[0] - Tension[0])/m;
ak[1] = (D[1] + L[1] - Tension[1] - m*g)/m;

vk[0] = vk[0] + h*ak[0];
vk[1] = vk[1] + h*ak[1];

rk[0] = rk[0] + h*vk[0];
rk[1] = rk[1] + h*vk[1];

// Check: Sum of total forces

Ftot[0] = L[0] + D[0] + Fg[0] + Tension[0];
Ftot[1] = L[1] + D[1] + Fg[1] + Tension[1];

#ifdef DEBUG
    printf("N=%f\n", N);
    printf("Vw[0]= %f, Vw[1]=%f\n", W[0], W[1]);
    printf("Vkx=%f, Vky=%f\n", vk[0], vk[1]); 
    printf("Va_mod=%f, va[0]=%f, va[1]=%fn", Va_mod, va[0], va[1]); 
    printf("L[0]=%f\n", L[0]); 
    printf("L[1]=%f\n", L[1]); 
    printf("D[0]=%f\n", D[0]); 
    printf("D[1]=%f\n", D[1]); 
    printf("Fg[0]=%f\n", Fg[0]); 
    printf("Fg[1]=%f\n", Fg[1]);
    printf("Tx=%f, Ty=%f\n", Tension[0], Tension[1]);
    printf("*theta = %f\n", *theta);
    printf("a_block[0] = %f, a_block[1] = %f\n", a_block[0], a_block[1]);
    printf("v_block[0] = %f, v_block[1] = %f\n", v_block[0], v_block[1]);
    printf("r_block[0] = %f, r_block[1] = %f\n", r_block[0], r_block[1]);
    printf("rk[0]= %f, rk[1]=%f\n", rk[0], rk[1]);
    printf("vk[0]= %f, vk[1]=%f\n", vk[0], vk[1]);
    printf("ak[0]= %f, ak[1]=%f\n", ak[0], ak[1]);
    printf("forze totali x:%f\n", Ftot[0]);
    printf("forze totali y:%f\n", Ftot[1]);
    printf("\n");
#endif    

}

#endif

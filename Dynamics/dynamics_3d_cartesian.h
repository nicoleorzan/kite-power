#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "constants.h"

#ifndef __dynamics__
#define __dynamics__

//#define DEBUG

#define h 0.00001                // Integration step
#define coeff_friction 0.4       // Block friction coefficient
#define m_block 10000.0          // Block mass

double W;
double Va_mod;
double T_denom;
double denom1, denom2;
double T1, T2;
double va[3];
double L[3], Lc;
double D[3], Dc;
double t1[3], t2[3], t3[3];
double t2_mod, t3_mod;
double Fg[3] = {0, 0, -m*g};
double F_aer[3];
double F_attrito[2];
double F_attrito_mod;
double N;
double Tension[3];
double r_diff_mod;

double T_orizz;
double F_attr_sign;
double Ftot[3]; 

void variables_initialization(double * rk, double * vk, double * ak, double theta, double phi,
                              double * r_block, double * v_block, double * a_block){
    r_block[0] = 0;
    r_block[1] = 0;
    r_block[2] = 0;

    v_block[0] = 0;
    v_block[1] = 0;
    v_block[2] = 0;

    a_block[0] = 0;
    a_block[1] = 0;
    a_block[2] = 0;

    rk[0] = r_block[0] + R*sin(theta)*cos(phi);
    rk[1] = r_block[1] + R*sin(theta)*sin(phi);
    rk[2] = R*cos(theta);

    vk[0] = 0;
    vk[1] = 0;
    vk[2] = 0;

    ak[0] = 0;
    ak[1] = 0;
    ak[2] = 0;
}

void integration_trajectory(double * rk, double * vk, double * ak, // Kite variables
                            double * r_block, double * v_block, double * a_block, // Block variables
                            double * r_diff, double * v_diff, double * a_diff, 
                            double * theta, double * phi,
                            int alpha, double mu,
                            double * W, double * lift, double * drag,
                            double * T, int it){

r_diff[0] = rk[0] - r_block[0];
r_diff[1] = rk[1] - r_block[1];
r_diff[2] = rk[2] - r_block[2];

//printf("rdiff0=%f, rdiff1=%f, rdiff2=%f\n", r_diff[0], r_diff[1], r_diff[2]);

*phi = atan(r_diff[1]/r_diff[0]);

*theta = acos(r_diff[2]/sqrt(r_diff[0]*r_diff[0] + r_diff[1]*r_diff[1] + r_diff[2]*r_diff[2] ));

r_diff_mod = sqrt(r_diff[0]*r_diff[0] + r_diff[1]*r_diff[1] + r_diff[2]*r_diff[2] );
                       
va[0] = vk[0] - W[0];              // Apparent velocity on x
va[1] = vk[1] - W[1];              // Apparent velocity on y
va[2] = vk[2] - W[2];              // Apparent velocity on z

Va_mod = sqrt(va[0]*va[0] + va[1]*va[1] + va[2]*va[2]);

// Compute vectorial products

t1[0] = sin(*theta)*cos(*phi);
t1[1] = sin(*theta)*sin(*phi);
t1[2] = cos(*theta);

// Computing t2

t2[0] = t1[1]*va[2] - t1[2]*va[1];
t2[1] = t1[2]*va[0] - t1[0]*va[2];
t2[2] = t1[0]*va[1] - t1[1]*va[0];

t2_mod = sqrt(t2[0]*t2[0] + t2[1]*t2[1] + t2[2]*t2[2]);

t2[0] = t2[0]/t2_mod;
t2[1] = t2[1]/t2_mod;
t2[2] = t2[2]/t2_mod;

// Computing t3

t3[0] = va[1]*t2[2] - va[2]*t2[1];
t3[1] = va[2]*t2[0] - va[0]*t2[2];
t3[2] = va[0]*t2[1] - va[1]*t2[0];

t3_mod = sqrt(t3[0]*t3[0] + t3[1]*t3[1] + t3[2]*t3[2]);

t3[0] = t3[0]/t3_mod;
t3[1] = t3[1]/t3_mod;
t3[2] = t3[2]/t3_mod;

// Computing Lift and Drag            

Lc = 0.5*rho*CL_alpha[alpha]*A*Va_mod*Va_mod;

L[0] = Lc*(t2[0]*sin(mu) + t3[0]*cos(mu));
L[1] = Lc*(t2[1]*sin(mu) + t3[1]*cos(mu));
L[2] = Lc*(t2[2]*sin(mu) + t3[2]*cos(mu));

Dc = 0.5*rho*CD_alpha[alpha]*A*Va_mod;

D[0] = -Dc*va[0];
D[1] = -Dc*va[1];
D[2] = -Dc*va[2];

F_aer[0] = L[0] + D[0];
F_aer[1] = L[1] + D[1];
F_aer[2] = L[2] + D[2];

// Solving motion equations:
// printf("Rdiffmod=%f\n", r_diff_mod);

denom1 = R*(m+m_block)/(m*m_block)
        - cos(*theta)/m_block*(rk[2] + coeff_friction*(cos(*phi)*r_diff[0] + sin(*phi)*r_diff[1]));

T1 = (F_aer[0]*r_diff[0] + F_aer[1]*r_diff[1] + F_aer[2]*r_diff[2])/m
    + (v_diff[0]*v_diff[0] + v_diff[1]*v_diff[1] + v_diff[2]*v_diff[2]) 
    - g*(r_diff[2] + coeff_friction*(cos(*phi)*r_diff[0] + sin(*phi)*r_diff[1]));

T1 = T1/denom1;

denom2 = R*(m+m_block)/(m*m_block)
    + cos(*theta)/m_block*(rk[2] + coeff_friction*(cos(*phi)*r_diff[0] + sin(*phi)*r_diff[1]));
T2 = (F_aer[0]*r_diff[0] + F_aer[1]*r_diff[1] + F_aer[2]*r_diff[2])/m
    + (v_diff[0]*v_diff[0] + v_diff[1]*v_diff[1] + v_diff[2]*v_diff[2]) 
    + g*(r_diff[2] + coeff_friction*(cos(*phi)*r_diff[0] + sin(*phi)*r_diff[1]));

T2 = T2/denom2;

*T = T1;

if (m_block*g < T2*cos(*theta)){
    printf("========>T2\n\n");
    //*T = T2;
}

/*if (m_block*g >= T1*cos(*theta)){
    *T = T1;
}
else if (m_block*g < T2*cos(*theta)){
    printf("========>T2\n\n");
    *T = T2;
}
else {
    printf("ERROR!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    printf("it=%d\n", it);
    printf("Mg=%f, T1cos=%f, T2cos=%f\n", m_block*g, T1*cos(*theta), T2*cos(*theta));
    printf("Mg-T1costheta = %f\n", m_block*g - T1*cos(*theta));
    printf("Mg-T2costheta = %f\n", m_block*g - T2*cos(*theta));
    //printf("T1=%f, T2=%f\n", T1, T2);
}*/

Tension[0] = *T*sin(*theta)*cos(*phi);
Tension[1] = *T*sin(*theta)*sin(*phi);
Tension[2] = *T*cos(*theta);

N = m_block*g - Tension[2];

// sulla forza d'attrito va messo il meno a mano?? forse si!
F_attrito_mod = coeff_friction*fabs(N);
F_attrito[0] = F_attrito_mod*cos(*phi);
F_attrito[1] = F_attrito_mod*sin(*phi);

a_block[0] = (Tension[0] + F_attrito[0])/m_block;
a_block[1] = (Tension[1] + F_attrito[1])/m_block;

v_block[0] = v_block[0] + h*a_block[0]; 
v_block[1] = v_block[1] + h*a_block[1];

r_block[0] = r_block[0] + h*v_block[0]; 
r_block[1] = r_block[1] + h*v_block[1];

ak[0] = (D[0] + L[0] - Tension[0])/m;
ak[1] = (D[1] + L[1] - Tension[1])/m;
ak[2] = (D[2] + L[2] - Tension[2] - m*g)/m;

vk[0] = vk[0] + h*ak[0];
vk[1] = vk[1] + h*ak[1];
vk[2] = vk[2] + h*ak[2];

rk[0] = rk[0] + h*vk[0];
rk[1] = rk[1] + h*vk[1];
rk[2] = rk[2] + h*vk[2];

// Check: Sum of total forces

Ftot[0] = L[0] + D[0] + Fg[0] + Tension[0];
Ftot[1] = L[1] + D[1] + Fg[1] + Tension[1];
Ftot[2] = L[2] + D[2] + Fg[2] + Tension[2];

#ifdef DEBUG
    printf("N=%f\n", N);
    printf("Vw[0]= %f, Vw[1]=%f, Vw[2]=%f\n", W[0], W[1], W[2]);
    printf("Vkx=%f, Vky=%f, Vkz=%f\n", vk[0], vk[1], vk[2]); 
    printf("Va_mod=%f, va[0]=%f, va[1]=%f, va[2]=%f\n", Va_mod, va[0], va[1], va[2]); 
    printf("L[0]=%f\n", L[0]); 
    printf("L[1]=%f\n", L[1]); 
    printf("L[2]=%f\n", L[2]); 
    printf("D[0]=%f\n", D[0]); 
    printf("D[1]=%f\n", D[1]); 
    printf("D[2]=%f\n", D[2]); 
    printf("Fg[0]=%f\n", Fg[0]); 
    printf("Fg[1]=%f\n", Fg[1]);
    printf("Fg[2]=%f\n", Fg[2]);
    printf("Tx=%f, Ty=%f, Tz=%f\n", Tension[0], Tension[1], Tension[2]);
    printf("*theta = %f\n", *theta);
    printf("*phi = %f\n", *phi);
    printf("a_block[0] = %f, a_block[1] = %f, a_block[2] = %f\n", a_block[0], a_block[1], a_block[2]);
    printf("v_block[0] = %f, v_block[1] = %f, v_block[2] = %f\n", v_block[0], v_block[1], v_block[2]);
    printf("r_block[0] = %f, r_block[1] = %f, r_block[2] = %f\n", r_block[0], r_block[1], r_block[2]);
    printf("rk[0]= %f, rk[1]=%f, rk[2]=%f\n", rk[0], rk[1], rk[2]);
    printf("vk[0]= %f, vk[1]=%f, vk[2]=%f\n", vk[0], vk[1], vk[2]);
    printf("ak[0]= %f, ak[1]=%f, ak[2]=%f\n", ak[0], ak[1], ak[2]);
    printf("forze totali x:%f\n", Ftot[0]);
    printf("forze totali y:%f\n", Ftot[1]);
    printf("forze totali z:%f\n", Ftot[2]);
    printf("\n");
#endif    

}

#endif

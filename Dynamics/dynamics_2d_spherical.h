#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "constants.h"

#ifndef __dynamics__
#define __dynamics__

//#define DEBUG

#define h 0.0001                // Integration step
#define coeff_friction 0.4      // Block friction coefficient
#define m_block 100.0          // Block mass
#define V_THRESHOLD 10E-6

double Va_mod;
double va[2];
double L[2];
double D[2];
double beta;
double Fg[2] = {0, -m*g};
double F_friction;
double N;
double Tension[2];
double Ftot[2]; 

double theta_star;

void variables_initialization(double * rk, double * vk, double * ak, 
                            double theta0, double vtheta0,
                            double * r_block, double * v_block, double * a_block,
                            double * theta){ // Angle, velocity and acceleration values
    r_block[0] = 0;
    r_block[1] = 0;

    v_block[0] = 0;
    v_block[1] = 0;

    a_block[0] = 0;
    a_block[1] = 0;

    theta[0] = theta0;
    theta[1] = vtheta0;
    theta[2] = 0.;

    rk[0] = R*cos(theta[0]) + r_block[0];
    rk[1] = R*sin(theta[0]);

    vk[0] = v_block[0] - R*sin(theta[0])*theta[1];
    vk[1] = R*cos(theta[0])*theta[1];

    ak[0] = a_block[0] - R*cos(theta[0])*theta[1]*theta[1] - R*sin(theta[0])*theta[2];
    ak[1] = -R*sin(theta[0])*theta[1]*theta[1] + R*cos(theta[0])*theta[2];
}

void integration_trajectory(double * rk, double * vk, double * ak, // Kite variables
                            double * r_block, double * v_block, double * a_block, // Block variables
                            double * theta, // Angle, velocity and acceleration values
                            double * T, int alpha,
                            double * W, double * lc, double * dc){

    va[0] = vk[0] - W[0];              // Apparent velocity on x
    va[1] = vk[1] - W[1];              // Apparent velocity on z

    Va_mod = sqrt(va[0]*va[0] + va[1]*va[1]);

    // Computing Lift and Drag

    beta = atan2(va[1], va[0]);

    *lc = 0.5*rho*CL_alpha[alpha]*A*Va_mod*Va_mod;

    *dc = 0.5*rho*CD_alpha[alpha]*A*Va_mod*Va_mod;

    // Lift (x, z)

    if (beta > - PI/2. && beta < PI/2.){
        L[0] = *lc*cos(beta + PI/2.);
        L[1] = *lc*sin(beta + PI/2.);
    } else {
        L[0] = *lc*cos(beta - PI/2.);
        L[1] = *lc*sin(beta - PI/2.);
    }

    // Drag (x, z)

    D[0] = *dc*cos(beta + PI);
    D[1] = *dc*sin(beta + PI);

    theta_star = atan((*lc - m*g)/(*dc));                 

    // Theta acceleration

    theta[2] = (a_block[0]*sin(theta[0]*m) - m*g*cos(theta[0]) - L[0]*sin(theta[0]) - D[0]*sin(theta[0]) + L[1]*cos(theta[0]) + D[1]*cos(theta[0]) )/(m*R);

    // Tension next step

    *T = - a_block[0]*cos(theta[0])*m + m*R*theta[1]*theta[1] - m*g*sin(theta[0]) + L[0]*cos(theta[0]) + L[1]*sin(theta[0]) + D[0]*cos(theta[0]) + D[1]*sin(theta[0]);

    Tension[0] = *T*cos(theta[0]);
    Tension[1] = *T*sin(theta[0]);

    // Compute Block motion

    N = m_block*g - Tension[1];

    F_friction = coeff_friction*fabs(N); // the sign is imposed into the equation below by hand

    if ( fabs(v_block[0]) < V_THRESHOLD ){ // blocco fermo
        if ( fabs(Tension[0]) > fabs(F_friction)  ){
            if (cos(theta[0]) > 0){
                a_block[0] = (Tension[0] - F_friction )/m_block;
            }
            else{
                a_block[0] = (Tension[0] + F_friction )/m_block;
            }
        }
        else { a_block[0] = 0; } 
    }
    else if ( v_block[0] > V_THRESHOLD ){  //v_block[0] > 0 ){  
        a_block[0] = (Tension[0] - N )/m_block;
    }
    else if ( v_block[0] < -V_THRESHOLD ){   //v_block[0] < 0 ){
        a_block[0] = (Tension[0] + N )/m_block;
    }

    a_block[1] = 0;

    v_block[0] = v_block[0] + h*a_block[0];
    v_block[1] = v_block[1] + h*a_block[1];

    r_block[0] = r_block[0] + h*v_block[0]; 
    r_block[1] = r_block[1] + h*v_block[1];

    theta[1] = theta[1] + h*theta[2];
    theta[0] = theta[0] + h*theta[1];

    // Kite position, velocity and acceleration update (x, z)

    ak[0] = a_block[0] - R*cos(theta[0])*theta[1]*theta[1] - R*sin(theta[0])*theta[2];
    ak[1] = -R*sin(theta[0])*theta[1]*theta[1] + R*cos(theta[0])*theta[2];

    vk[0] = v_block[0] - R*sin(theta[0])*theta[1];
    vk[1] = R*cos(theta[0])*theta[1];
    
    rk[0] = r_block[0] + R*cos(theta[0]);
    rk[1] = R*sin(theta[0]);

    // Check: Sum of total forces on kite

    Ftot[0] = L[0] + D[0] + Fg[0] - Tension[0];

    Ftot[1] = L[1] + D[1] + Fg[1] - Tension[1];

    // Power

    //*power = *T*dpos[0]*sin(pos[1]);

    #ifdef DEBUG
        printf("thetastar = %f\n", theta_star);
        printf("beta=%f\n", beta);
        printf("W[0]= %f, W[1]=%f\n", W[0], W[1]);
        printf("Vkx=%f, Vkz=%f\n", vk[0], vk[1]); 
        printf("Va_mod=%f, va[0]=%f, va[1]=%f\n", Va_mod, va[0], va[1]); 
        printf("L[0]=%f\n", L[0]); 
        printf("L[1]=%f\n", L[1]); 
        printf("D[0]=%f\n", D[0]); 
        printf("D[1]=%f\n", D[1]); 
        printf("Fg[0]=%f\n", Fg[0]); 
        printf("Fg[1]=%f\n", Fg[1]); 
        printf("T=%f, a_block[0]=%f, acc_theta = %f\n", *T, a_block[0], theta[2]);
        printf("Tx=%f, Ty=%f\n", *T*cos(theta[0]), *T*sin(theta[0]));
        printf("theta[0] = %f, theta[1] = %f\n", theta[0], theta[1]);
        printf("a_block[0] = %f, a_block[1] = %f\n", a_block[0], a_block[1]);
        printf("v_block[0] = %f, v_block[1] = %f\n", v_block[0], v_block[1]);
        printf("r_block[0] = %f, r_block[1] = %f\n", r_block[0], r_block[1]);
        printf("rk[0]= %f, rk[1]=%f\n", rk[0], rk[1]);
        printf("vk[0]= %f, vk[1]=%f\n", vk[0], vk[1]);
        printf("ak[0]= %f, ak[1]=%f\n", ak[0], ak[1]);
        printf("forze totali x:%f\n", Ftot[0]);
        printf("forze totali x:%f\n", Ftot[1]);
        printf("\n");
    #endif    

}

#endif

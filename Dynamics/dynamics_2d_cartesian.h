#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "constants.h"

#ifndef __dynamics__
#define __dynamics__

//#define DEBUG

#define coeff_friction 0.4      // Block friction coefficient
#define V_THRESHOLD 10E-6

double Va_mod;
double va[2];
double L[2];
double D[2];
double beta;
double Fg[2] = {0, -m*g};
double F_friction;
double N;
double Tension[2], denom;
double Ftot[2]; 

double F_aer[2];

void variables_initialization(double * rk, double * vk, double * ak,
                             double theta, double dtheta,
                             double * r_block, double * v_block, double * a_block){
    r_block[0] = 0;
    r_block[1] = 0;

    v_block[0] = 0;
    v_block[1] = 0;

    a_block[0] = 0;
    a_block[1] = 0;

    rk[0] = r_block[0] + R*cos(theta);
    rk[1] = r_block[1] + R*sin(theta);

    vk[0] = -R*dtheta*sin(theta);
    vk[1] = R*dtheta*cos(theta);

    ak[0] = -R*dtheta*dtheta*cos(theta);
    ak[1] = -R*dtheta*dtheta*sin(theta);
}

void integration_trajectory(double * rk, double * vk, double * ak, // Kite variables
                            double * r_block, double * v_block, double * a_block, // Block variables
                            double * r_diff, double * v_diff, double * a_diff, 
                            double * theta,
                            int alpha,
                            double * W, double * lc, double * dc,
                            double * T, int it){

    r_diff[0] = rk[0] - r_block[0];
    r_diff[1] = rk[1] - r_block[1];

    v_diff[0] = vk[0] - v_block[0];
    v_diff[1] = vk[1] - v_block[1];

    a_diff[0] = ak[0] - a_block[0];
    a_diff[1] = ak[1] - a_block[1];

    *theta = atan2(r_diff[1], r_diff[0]);
                        
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

    F_aer[0] = L[0] + D[0];
    F_aer[1] = L[1] + D[1];

    // Computing block motion

    if ( fabs(v_block[0]) < V_THRESHOLD ){ // block not moving
        
        denom = R*(m+m_block)/(m*m_block)
                - sin(*theta)/m_block*(r_diff[1] - coeff_friction*cos(*theta)*r_diff[0]);

        *T = (F_aer[0]*r_diff[0] + F_aer[1]*r_diff[1])/m
            + (v_diff[0]*v_diff[0] + v_diff[1]*v_diff[1]) 
            - g*(r_diff[1] - coeff_friction*cos(*theta)*r_diff[0]);

        *T = *T/denom;

        Tension[0] = *T*cos(*theta);
        Tension[1] = *T*sin(*theta);
    
        N = m_block*g - Tension[1];

        F_friction = -coeff_friction*fabs(N)*cos(*theta);

        #ifdef DEBUG
            printf("v almost zero\nv_block=%.10f\n", v_block[0]);
            printf("*T=%f, theta=%f\n", *T, *theta);
            printf("Tx=%f, Ty=%f\n", Tension[0], Tension[1]);
            printf("F_friction=%f\n", F_friction);
            printf("\n");
        #endif

        if ( fabs(Tension[0]) > fabs(F_friction) ){
            a_block[0] = ( Tension[0] + F_friction )/m_block;
        }
        else { // F_friction = -Tension[0], ==> a_block[0] = 0
            
            denom = R*(m+m_block)/(m*m_block) 
            - sin(*theta)*r_diff[1]/m_block - cos(*theta)*r_diff[0]/m_block;

            *T = (F_aer[0]*r_diff[0] + F_aer[1]*r_diff[1])/m 
                + (v_diff[0]*v_diff[0] + v_diff[1]*v_diff[1]) - g*r_diff[1];

            Tension[0] = *T*cos(*theta);
            Tension[1] = *T*sin(*theta);

            F_friction = -Tension[0];
            
            a_block[0] = ( Tension[0] + F_friction )/m_block;
        } 
    }
    else { // block moving, |v| > 10E-6

        denom = R*(m+m_block)/(m*m_block)
                - sin(*theta)/m_block*(r_diff[1] - coeff_friction*r_diff[0]*v_block[0]/fabs(v_block[0]));

        *T = (F_aer[0]*r_diff[0] + F_aer[1]*r_diff[1])/m
            + (v_diff[0]*v_diff[0] + v_diff[1]*v_diff[1]) 
            - g*(r_diff[1] - coeff_friction*r_diff[0]*v_block[0]/fabs(v_block[0]));

        *T = *T/denom;

        Tension[0] = *T*cos(*theta);
        Tension[1] = *T*sin(*theta);
    
        N = m_block*g - Tension[1];

        F_friction = -coeff_friction*fabs(N)*v_block[0]/fabs(v_block[0]);

        a_block[0] = ( Tension[0] + F_friction )/m_block;

    #ifdef DEBUG
        printf("PARTE2\nv_block=%.10f\n", v_block[0]);
        printf("*T=%f, theta=%f\n", *T, *theta);
        printf("Tx=%f, Ty=%f\n", Tension[0], Tension[1]);
        printf("F_friction=%f\n", F_friction);
        printf("\n");
    #endif
    }
    //printf("theta=%f, *T=%f\n", *theta,*T);

    a_block[1] = 0;

    v_block[0] = v_block[0] + h*a_block[0]; 
    v_block[1] = v_block[1] + h*a_block[1];

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
        printf("Va_mod=%f, va[0]=%f, va[1]=%f\n", Va_mod, va[0], va[1]); 
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

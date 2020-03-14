#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "constants.h"

#ifndef __dynamics__
#define __dynamics__

//#define DEBUG

#define coeff_friction 0.4       // Block friction coefficient
#define V_THRESHOLD 10E-6

double Va_mod;
double va[3];
double L[3];
double D[3];
double F_aer[3];
double Fg[3] = {0, 0, -m*g};
double F_friction[2];
double N;
double Tension[3], denom;
double Ftot[3];

double t1[3], t2[3], t3[3];
double t2_mod, t3_mod;
double v_block_mod;

void variables_initialization(double * rk, double * vk, double * ak, double theta, double phi,
                             double vtheta, double vphi, double * r_block, double * v_block, double * a_block){
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

    vk[0] = v_block[0] + R*cos(theta)*cos(phi)*vtheta - R*sin(theta)*sin(phi)*vphi;
    vk[1] = v_block[1] + R*cos(theta)*sin(phi)*vtheta + R*sin(theta)*cos(phi)*vphi;
    vk[2] = -R*sin(theta)*vtheta;

    ak[0] = R*(a_block[0]/R -sin(theta)*cos(phi)*(vtheta*vtheta + vphi*vphi)
            -2*cos(theta)*sin(phi)*vtheta*vphi);
    ak[1] = R*(a_block[1]/R - sin(theta)*sin(phi)*(vtheta*vtheta + vphi*vphi)
            +2*cos(theta)*cos(phi)*vtheta*vphi);
    ak[2] = -R*cos(theta)*vtheta*vtheta;
}

void integration_trajectory(double * rk, double * vk, double * ak, // Kite variables
                            double * r_block, double * v_block, double * a_block, // Block variables
                            double * r_diff, double * v_diff, double * a_diff, 
                            double * theta, double * phi,
                            int alpha, double mu,
                            double * W, double * lc, double * dc,
                            double * T, int it){

    r_diff[0] = rk[0] - r_block[0];
    r_diff[1] = rk[1] - r_block[1];
    r_diff[2] = rk[2] - r_block[2];

    v_diff[0] = vk[0] - v_block[0];
    v_diff[1] = vk[1] - v_block[1];
    v_diff[2] = vk[2] - v_block[2];

    a_diff[0] = ak[0] - a_block[0];
    a_diff[1] = ak[1] - a_block[1];
    a_diff[2] = ak[2] - a_block[2];

    printf("theta=%f, ", *theta);
    printf("phi=%f\n", *phi);

    *phi = atan(r_diff[1]/r_diff[0]); // CHECK THE ANGLES!!!!!!!!!!!!!!!!!!!!!!!! <=========

    *theta = acos(r_diff[2]/sqrt(r_diff[0]*r_diff[0] + r_diff[1]*r_diff[1] + r_diff[2]*r_diff[2] ));
                        
    va[0] = vk[0] - W[0];              // Apparent velocity on x
    va[1] = vk[1] - W[1];              // Apparent velocity on y
    va[2] = vk[2] - W[2];              // Apparent velocity on z

    Va_mod = sqrt(va[0]*va[0] + va[1]*va[1] + va[2]*va[2]);

    // Computing Lift and Drag   

    // computing t1 

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

    *lc = 0.5*rho*CL_alpha[alpha]*A*Va_mod*Va_mod;

    L[0] = *lc*(t2[0]*sin(mu) + t3[0]*cos(mu));
    L[1] = *lc*(t2[1]*sin(mu) + t3[1]*cos(mu));
    L[2] = *lc*(t2[2]*sin(mu) + t3[2]*cos(mu));

    *dc = 0.5*rho*CD_alpha[alpha]*A*Va_mod;

    D[0] = *dc*va[0];
    D[1] = *dc*va[1];
    D[2] = *dc*va[2];

    F_aer[0] = L[0] + D[0];
    F_aer[1] = L[1] + D[1];
    F_aer[2] = L[2] + D[2];

    // Solving dynamics

    if ( sqrt(v_block[0]*v_block[0] + v_block[1]*v_block[1]) < V_THRESHOLD ){ // block not moving

        denom = R*(m+m_block)/(m*m_block)
            + cos(*theta)/m_block*(coeff_friction*(sin(*theta)*cos(*phi)*r_diff[0] + sin(*theta)*sin(*phi)*r_diff[1]) - r_diff[2]);

        *T = (F_aer[0]*r_diff[0] + F_aer[1]*r_diff[1] + F_aer[2]*r_diff[2])/m
        + (v_diff[0]*v_diff[0] + v_diff[1]*v_diff[1] + v_diff[2]*v_diff[2])
        + g*(coeff_friction/m_block*(sin(*theta)*cos(*phi)*r_diff[0] + sin(*theta)*sin(*phi)*r_diff[1]) -r_diff[2]);

        *T = *T/denom;

        Tension[0] = *T*sin(*theta)*cos(*phi);
        Tension[1] = *T*sin(*theta)*sin(*phi);
        Tension[2] = *T*cos(*theta);

        N = m_block*g - Tension[2];

        F_friction[0] = -coeff_friction*fabs(N)*sin(*theta)*cos(*phi);
        F_friction[1] = -coeff_friction*fabs(N)*sin(*theta)*sin(*phi);

        if ( fabs(Tension[0]) >= fabs(F_friction[0]) && fabs(Tension[1]) >= fabs(F_friction[1]) ){ // both bigger

            a_block[0] = ( Tension[0] + F_friction[0] )/m_block;
            a_block[1] = ( Tension[1] + F_friction[1] )/m_block;

        } else if (  fabs(Tension[0]) < fabs(F_friction[0]) && fabs(Tension[1]) < fabs(F_friction[1])  ) { // both smaller

            denom = R*(m+m_block)/(m*m_block)
            - 1/m_block*( sin(*theta)*( cos(*phi)*r_diff[0] + sin(*phi)*r_diff[1] ) - cos(*theta)*r_diff[2] );

            *T = (F_aer[0]*r_diff[0] + F_aer[1]*r_diff[1] + F_aer[2]*r_diff[2])/m
            + (v_diff[0]*v_diff[0] + v_diff[1]*v_diff[1] + v_diff[2]*v_diff[2])
            -g*r_diff[2];

            *T = *T/denom;

            Tension[0] = *T*sin(*theta)*cos(*phi);
            Tension[1] = *T*sin(*theta)*sin(*phi);
            Tension[2] = *T*cos(*theta);

            N = m_block*g - Tension[2];

            F_friction[0] = -Tension[0];
            F_friction[1] = -Tension[1];

            a_block[0] = ( Tension[0] + F_friction[0] )/m_block;
            a_block[1] = ( Tension[1] + F_friction[1] )/m_block;

        } else if (  fabs(Tension[0]) >= fabs(F_friction[0]) && fabs(Tension[1]) < fabs(F_friction[1])  ) { 

            denom = R*(m+m_block)/(m*m_block)
            + 1/m_block*(coeff_friction*cos(*theta)*sin(*theta)*cos(*phi)*r_diff[0] 
            - sin(*theta)*sin(*phi)*r_diff[1] - cos(*theta)*r_diff[2]);

            *T = (F_aer[0]*r_diff[0] + F_aer[1]*r_diff[1] + F_aer[2]*r_diff[2])/m
            + (v_diff[0]*v_diff[0] + v_diff[1]*v_diff[1] + v_diff[2]*v_diff[2])
            + coeff_friction*g*sin(*theta)*cos(*phi)*r_diff[0];

            *T = *T/denom;

            Tension[0] = *T*sin(*theta)*cos(*phi);
            Tension[1] = *T*sin(*theta)*sin(*phi);
            Tension[2] = *T*cos(*theta);

            N = m_block*g - Tension[2];

            F_friction[0] = -coeff_friction*fabs(N)*sin(*theta)*cos(*phi);
            F_friction[1] = -Tension[1];

            a_block[0] = ( Tension[0] + F_friction[0] )/m_block;
            a_block[1] = ( Tension[1] + F_friction[1] )/m_block;

        } else if (  fabs(Tension[0]) < fabs(F_friction[0]) && fabs(Tension[1]) >= fabs(F_friction[1])  ) { 

            denom = R*(m+m_block)/(m*m_block)
            + 1/m_block*(coeff_friction*cos(*theta)*sin(*theta)*sin(*phi)*r_diff[1] 
            - sin(*theta)*cos(*phi)*r_diff[0] - cos(*theta)*r_diff[2]);

            *T = (F_aer[0]*r_diff[0] + F_aer[1]*r_diff[1] + F_aer[2]*r_diff[2])/m
            + (v_diff[0]*v_diff[0] + v_diff[1]*v_diff[1] + v_diff[2]*v_diff[2])
            + coeff_friction*g*sin(*theta)*sin(*phi)*r_diff[1];

            *T = *T/denom;

            Tension[0] = *T*sin(*theta)*cos(*phi);
            Tension[1] = *T*sin(*theta)*sin(*phi);
            Tension[2] = *T*cos(*theta);

            N = m_block*g - Tension[2];

            F_friction[0] = -Tension[0];
            F_friction[1] = -coeff_friction*fabs(N)*sin(*theta)*sin(*phi);

            a_block[0] = ( Tension[0] + F_friction[0] )/m_block;
            a_block[1] = ( Tension[1] + F_friction[1] )/m_block;
        }

        //printf("v < threshold: v=%f, T=%f, Tx=%f, Ty=%f\n", v_block_mod, *T,  Tension[0], Tension[1]);

    }

    else { // block moving, |v| > 10E-6

        v_block_mod = sqrt(v_block[0]*v_block[0] + v_block[1]*v_block[1]);

        denom = R*(m+m_block)/(m*m_block) + cos(*theta)/m_block*(-r_diff[2] + 
        coeff_friction/v_block_mod*(v_block[0]*r_diff[0] + v_block[1]*r_diff[1]) );

        *T = (F_aer[0]*r_diff[0] + F_aer[1]*r_diff[1] + F_aer[2]*r_diff[2])/m
        + (v_diff[0]*v_diff[0] + v_diff[1]*v_diff[1] + v_diff[2]*v_diff[2])
        + g*(-r_diff[2] + coeff_friction/v_block_mod*(v_block[0]*r_diff[0] + v_block[1]*r_diff[1]) );

        *T = *T/denom;

        Tension[0] = *T*sin(*theta)*cos(*phi);
        Tension[1] = *T*sin(*theta)*sin(*phi);
        Tension[2] = *T*cos(*theta);

        N = m_block*g - Tension[2];

        F_friction[0] = -coeff_friction*fabs(N)*v_block[0]/v_block_mod;
        F_friction[1] = -coeff_friction*fabs(N)*v_block[1]/v_block_mod;

        a_block[0] = ( Tension[0] + F_friction[0] )/m_block;
        a_block[1] = ( Tension[1] + F_friction[1] )/m_block;

        //printf("v > threshold: v=%f, T=%f, Tx=%f, Ty=%f\n", v_block_mod, *T,  Tension[0], Tension[1]);
    }

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

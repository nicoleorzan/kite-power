#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "constants1.h"

#ifndef __dynamics__
#define __dynamics__

#define CONST 10000

//#define DEBUG

double Va_mod;
double va[3];
double L[3];
double D[3];
double F_aer[3];
double Fg[3] = {0, 0, -m*g};
double F_friction[2];
double N;
double T1, denom1;
double T2, denom2;
double Tension[3], denom;
double Ftot[3];

double t1[3], t2[3], t3[3];
double t2_mod, t3_mod;
double v_block_mod;

void variables_initialization(double *rk, double *vk, double *ak, 
                             double theta, double phi,
                             double dtheta, double dphi, 
                             double *r_block, double *v_block, double *a_block,
                             double *r_diff, double *v_diff, double *a_diff){
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

    vk[0] = v_block[0] + R*cos(theta)*cos(phi)*dtheta - R*sin(theta)*sin(phi)*dphi;
    vk[1] = v_block[1] + R*cos(theta)*sin(phi)*dtheta + R*sin(theta)*cos(phi)*dphi;
    vk[2] = -R*sin(theta)*dtheta;
    //printf("vk[0]=%f, vk[1]=%f, vk[2]=%f\n", vk[0], vk[1], vk[2]);

    ak[0] = R*(a_block[0]/R -sin(theta)*cos(phi)*(dtheta*dtheta + dphi*dphi)
            -2*cos(theta)*sin(phi)*dtheta*dphi);
    ak[1] = R*(a_block[1]/R - sin(theta)*sin(phi)*(dtheta*dtheta + dphi*dphi)
            +2*cos(theta)*cos(phi)*dtheta*dphi);
    ak[2] = -R*cos(theta)*dtheta*dtheta;

    r_diff[0] = rk[0] - r_block[0];
    r_diff[1] = rk[1] - r_block[1];
    r_diff[2] = rk[2] - r_block[2];

    v_diff[0] = vk[0] - v_block[0];
    v_diff[1] = vk[1] - v_block[1];
    v_diff[2] = vk[2] - v_block[2];

    a_diff[0] = ak[0] - a_block[0];
    a_diff[1] = ak[1] - a_block[1];
    a_diff[2] = ak[2] - a_block[2];
}

void integration_trajectory(double * rk, double * vk, double * ak, // Kite variables
                            double * r_block, double * v_block, double * a_block, // Block variables
                            double * r_diff, double * v_diff, double * a_diff, 
                            double * theta, double * phi,
                            int alpha, double mu,
                            double * W, double * lc, double * dc,
                            double * T, double *F_attr, int it, int * sector, 
                            double * l0, double * l1, double * l2, double * d0, double * d1, double * d2){

    //printf("\ni=%d, noi\n", it);                               
    
    r_diff[0] = rk[0] - r_block[0];
    r_diff[1] = rk[1] - r_block[1];
    r_diff[2] = rk[2] - r_block[2];

    v_diff[0] = vk[0] - v_block[0];
    v_diff[1] = vk[1] - v_block[1];
    v_diff[2] = vk[2] - v_block[2];

    a_diff[0] = ak[0] - a_block[0];
    a_diff[1] = ak[1] - a_block[1];
    a_diff[2] = ak[2] - a_block[2];

    *theta = atan2(sqrt(r_diff[0]*r_diff[0] + r_diff[1]*r_diff[1]), r_diff[2]);
    *phi = atan2(r_diff[1], r_diff[0]);
                       
    va[0] = vk[0] - W[0];              // Apparent velocity on x
    va[1] = vk[1] - W[1];              // Apparent velocity on y
    va[2] = vk[2] - W[2];              // Apparent velocity on z

    Va_mod = sqrt(va[0]*va[0] + va[1]*va[1] + va[2]*va[2]);
    //printf("va0=%f, va1=%f, va2=%f\n", va[0]/Va_mod, va[1]/Va_mod, va[2]/Va_mod);

    // Computing Lift and Drag   

    // computing t1 
    
    t1[0] = sin(*theta)*cos(*phi);
    t1[1] = sin(*theta)*sin(*phi);
    t1[2] = cos(*theta);
    //printf("t1[0]=%f, t1[1]=%f, t1[2]=%f\n", t1[0], t1[1], t1[2]);

    t2[0] = (t1[1]*va[2] - t1[2]*va[1]);
    t2[1] = (t1[2]*va[0] - t1[0]*va[2]);
    t2[2] = (t1[0]*va[1] - t1[1]*va[0]); 

    t2_mod = sqrt(t2[0]*t2[0] + t2[1]*t2[1] + t2[2]*t2[2]);

    t2[0] = t2[0]/t2_mod;
    t2[1] = t2[1]/t2_mod;
    t2[2] = t2[2]/t2_mod;
    //printf("t20=%f, t21=%f, t22=%f\n", t2[0], t2[1], t2[2]);

    // Computing t3 = vrel x t2

    t3[0] = va[1]*t2[2] - va[2]*t2[1];
    t3[1] = va[2]*t2[0] - va[0]*t2[2];
    t3[2] = va[0]*t2[1] - va[1]*t2[0];

    t3_mod = sqrt(t3[0]*t3[0] + t3[1]*t3[1] + t3[2]*t3[2]);

    t3[0] = t3[0]/t3_mod;
    t3[1] = t3[1]/t3_mod;
    t3[2] = t3[2]/t3_mod;
    //printf("t30=%f, t31=%f, t32=%f\n", t3[0], t3[1], t3[2]);
    
    // Computing Lift and Drag          

    *lc = 0.5*rho*CL_alpha[alpha]*A*Va_mod*Va_mod;
    L[0] = *lc*(t2[0]*sin(mu) + t3[0]*cos(mu));
    L[1] = *lc*(t2[1]*sin(mu) + t3[1]*cos(mu));
    L[2] = *lc*(t2[2]*sin(mu) + t3[2]*cos(mu));
    //printf("L[0]=%f, L[1]=%f, L[2]=%f\n", L[0], L[1], L[2]);

    *l0 = L[0];
    *l1 = L[1];
    *l2 = L[2];

    *dc = 0.5*rho*CD_alpha[alpha]*A*Va_mod;

    D[0] = -*dc*va[0];
    D[1] = -*dc*va[1];
    D[2] = -*dc*va[2];
    //printf("D[0]=%f, D[1]=%f, D[2]=%f\n\n", D[0], D[1], D[2]);

    *d0 = D[0];
    *d1 = D[1];
    *d2 = D[2];

    *dc = 0.5*rho*CD_alpha[alpha]*A*Va_mod*sqrt(va[0]*va[0] + va[1]*va[1] + va[2]*va[2]);

    F_aer[0] = L[0] + D[0];
    F_aer[1] = L[1] + D[1];
    F_aer[2] = L[2] + D[2];

    // Compute tension

    // ===================== CASE 1) BLOCK NOT MOVING ( |v-block| < 10E-6 ) ==================

    v_block_mod = sqrt(v_block[0]*v_block[0] + v_block[1]*v_block[1]);

    if ( v_block_mod < V_THRESHOLD ){

        // |Mg| > |Tz|

        denom1 = R*(m+m_block)/(m*m_block)
            - cos(*theta)/m_block*( r_diff[2] - coeff_friction*(cos(*phi)*r_diff[0] + sin(*phi)*r_diff[1]));

        T1 = (F_aer[0]*r_diff[0] + F_aer[1]*r_diff[1] + F_aer[2]*r_diff[2])/m
        + (v_diff[0]*v_diff[0] + v_diff[1]*v_diff[1] + v_diff[2]*v_diff[2])
        - g*(r_diff[2] - coeff_friction*(cos(*phi)*r_diff[0] + sin(*phi)*r_diff[1]));

        T1 = T1/denom1;

        // |Mg| < |Tz|

        denom2 = R*(m+m_block)/(m*m_block)
            - cos(*theta)/m_block*( r_diff[2] + coeff_friction*(cos(*phi)*r_diff[0] + sin(*phi)*r_diff[1]));

        T2 = (F_aer[0]*r_diff[0] + F_aer[1]*r_diff[1] + F_aer[2]*r_diff[2])/m
        + (v_diff[0]*v_diff[0] + v_diff[1]*v_diff[1] + v_diff[2]*v_diff[2])
        - g*(r_diff[2] + coeff_friction*(cos(*phi)*r_diff[0] + sin(*phi)*r_diff[1]));

        T2 = T2/denom2;

        if ( m_block*g > T1*cos(*theta) ){
            *sector = 1;
            *T = T1;
        } else if ( m_block*g <= T2*cos(*theta) ){
            *sector = 2;
            *T = T2;
        } else { 
            printf("ERROR!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"); 
        }

        Tension[0] = *T*sin(*theta)*cos(*phi);
        Tension[1] = *T*sin(*theta)*sin(*phi);
        Tension[2] = *T*cos(*theta);

        N = m_block*g - Tension[2];

        F_friction[0] = -coeff_friction*fabs(N)*cos(*phi);
        F_friction[1] = -coeff_friction*fabs(N)*sin(*phi);

        *F_attr = sqrt(F_friction[0]*F_friction[0] + F_friction[1]*F_friction[1]);

        a_block[0] = (Tension[0] + F_friction[0])/m_block;
        a_block[1] = (Tension[1] + F_friction[1])/m_block;

        // ============> If the computed tension is bigger than friction force, block moves (so we have to do nothing).

        // ============> If not, recompute friction as: F_friction = -Tension[0]; which gives a_block[0] = 0

        if ( sqrt(Tension[0]*Tension[0] + Tension[1]*Tension[1]) < sqrt(F_friction[0]*F_friction[0] + F_friction[1]*F_friction[1]) ){

            *sector = 3;
            
            denom = R*(m+m_block)/(m*m_block) 
            - sin(*theta)/m_block*(cos(*phi)*r_diff[0] + sin(*phi)*r_diff[1]) - cos(*theta)/m_block*r_diff[2];

            *T = (F_aer[0]*r_diff[0] + F_aer[1]*r_diff[1] + F_aer[2]*r_diff[2])/m
            + (v_diff[0]*v_diff[0] + v_diff[1]*v_diff[1] + v_diff[2]*v_diff[2]) - g*r_diff[2];

            *T = *T/denom;

            Tension[0] = *T*sin(*theta)*cos(*phi);
            Tension[1] = *T*sin(*theta)*sin(*phi);
            Tension[2] = *T*cos(*theta);

            F_friction[0] = -Tension[0];
            F_friction[1] = -Tension[1];

            *F_attr = sqrt(F_friction[0]*F_friction[0] + F_friction[1]*F_friction[1]);

            a_block[0] = ( Tension[0] + F_friction[0] )/m_block;
            a_block[0] = ( Tension[1] + F_friction[1] )/m_block;

        }
    }

    // ========================== CASE 2) BLOCK MOVING (|v| >= 10E-6) ===> Fmu = -mu*|N|*vx/|vx| =====================

    else { 
        
        // |Mg| > |Tz|

        v_block_mod = sqrt(v_block[0]*v_block[0] + v_block[1]*v_block[1]);

        denom1 = R*(m+m_block)/(m*m_block) - cos(*theta)/m_block*(r_diff[2] - 
            coeff_friction/v_block_mod*(v_block[0]*r_diff[0] + v_block[1]*r_diff[1]) );

        T1 = (F_aer[0]*r_diff[0] + F_aer[1]*r_diff[1] + F_aer[2]*r_diff[2])/m
            + (v_diff[0]*v_diff[0] + v_diff[1]*v_diff[1] + v_diff[2]*v_diff[2])
            - g*(r_diff[2] - coeff_friction/v_block_mod*(v_block[0]*r_diff[0] + v_block[1]*r_diff[1]) );

        T1 = T1/denom1;

        // |Mg| < |Tz|

        denom2 = R*(m+m_block)/(m*m_block) - cos(*theta)/m_block*(r_diff[2] + 
            coeff_friction/v_block_mod*(v_block[0]*r_diff[0] + v_block[1]*r_diff[1]) );

        T2 = (F_aer[0]*r_diff[0] + F_aer[1]*r_diff[1] + F_aer[2]*r_diff[2])/m
            + (v_diff[0]*v_diff[0] + v_diff[1]*v_diff[1] + v_diff[2]*v_diff[2])
            - g*(r_diff[2] + coeff_friction/v_block_mod*(v_block[0]*r_diff[0] + v_block[1]*r_diff[1]) );

        T2 = T2/denom2;

        if ( m_block*g > T1*cos(*theta) ) {
            *sector = 4;
            *T = T1;
        } else if ( m_block*g <= T2*cos(*theta) ){
            *sector = 5;
            *T = T2;
        } else {
            printf("ERROR!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"); 
        }

        Tension[0] = *T*sin(*theta)*cos(*phi);
        Tension[1] = *T*sin(*theta)*sin(*phi);
        Tension[2] = *T*cos(*theta);

        N = m_block*g - Tension[2];

        F_friction[0] = -coeff_friction*fabs(N)*v_block[0]/v_block_mod;
        F_friction[1] = -coeff_friction*fabs(N)*v_block[1]/v_block_mod;

        *F_attr = sqrt(F_friction[0]*F_friction[0] + F_friction[1]*F_friction[1]);

        a_block[0] = ( Tension[0] + F_friction[0] )/m_block;
        a_block[1] = ( Tension[1] + F_friction[1] )/m_block;

    }

    a_block[2] = 0;

    // BLOCKED MOTION ON Y BECAUSE OF THE RAIL 
    
    v_block[0] = v_block[0] + h*a_block[0]; 
    v_block[1] = 0;//v_block[1] + h*a_block[1];

    r_block[0] = r_block[0] + h*v_block[0]; 
    r_block[1] = 0;//r_block[1] + h*v_block[1];

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

}

#endif

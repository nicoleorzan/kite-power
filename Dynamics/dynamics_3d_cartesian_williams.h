#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "constants.h"

#ifndef __dynamics__
#define __dynamics__

#define CONST 10000

//#define DEBUG

double Va_mod;
double va[3];
double va1[3];
double L[3];
double fabsL1[3];
double D[3];
double F_aer[3];
double Fg[3] = {0, 0, -m*g};
double F_friction[2];
double N;
double T1, denom1;
double T2, denom2;
double Tension[3], denom;
double Ftot[3];
double Txy_norm;
double L_mod;

double t1[3], t2[3], t3[3];
double t2_mod, t3_mod;
double v_block_mod;

double r_diff_modulo, v_diff_modulo;

double vk_spher[3];
double scal_norm;
double mod_sqrt;
double wep_mod;
double e0[3];
double wep[3];
double ew[3];
double W_spher[3];
double we[3], ewep[3], et[3];
double nwep, nwe, eta;
double Caer;
double Psi = 0;
double F_aerh[3];
double dtheta, dphi;
double L_mod, L_mod1;
double L1[3], L1c[3];
double D1[3], D1c[3];

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
    rk[1] = r_block[1] + R*sin(phi);
    rk[2] = R*cos(theta)*cos(phi);

    vk[0] = v_block[0] - R*sin(theta)*sin(phi)*dphi + R*cos(theta)*cos(phi)*dtheta;
    vk[1] = v_block[1] + R*cos(phi)*dphi;
    vk[2] = -R*cos(phi)*cos(theta)*dphi - cos(phi)*cos(theta)*dtheta;
    //printf("vk[0]=%f, vk[1]=%f, vk[2]=%f\n", vk[0], vk[1], vk[2]);

    ak[0] = 0; // da rifareR*(a_block[0]/R -sin(theta)*cos(phi)*(dtheta*dtheta + dphi*dphi)
            //-2*cos(theta)*sin(phi)*dtheta*dphi);
    ak[1] = 0;//R*(a_block[1]/R - sin(theta)*sin(phi)*(dtheta*dtheta + dphi*dphi)
            //+2*cos(theta)*cos(phi)*dtheta*dphi);
    ak[2] = 0;//-R*cos(theta)*dtheta*dtheta;

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
                            double *l0, double *l1, double *l2, double *d0, double *d1, double *d2){

    r_diff[0] = rk[0] - r_block[0];
    r_diff[1] = rk[1] - r_block[1];
    r_diff[2] = rk[2] - r_block[2];

    v_diff[0] = vk[0] - v_block[0];
    v_diff[1] = vk[1] - v_block[1];
    v_diff[2] = vk[2] - v_block[2];

    a_diff[0] = ak[0] - a_block[0];
    a_diff[1] = ak[1] - a_block[1];
    a_diff[2] = ak[2] - a_block[2];

    *phi = atan2(r_diff[1], sqrt(r_diff[0]*r_diff[0] + r_diff[1]*r_diff[1])); //atan2(r_diff[1], r_diff[0]);

    *theta = atan2(r_diff[0], r_diff[2]); //atan2(sqrt(r_diff[0]*r_diff[0] + r_diff[1]*r_diff[1]), r_diff[2]);
                        
    va[0] = vk[0] - W[0];              // Apparent velocity on x
    va[1] = vk[1] - W[1];              // Apparent velocity on y
    va[2] = vk[2] - W[2];              // Apparent velocity on z

    Va_mod = sqrt(va[0]*va[0] + va[1]*va[1] + va[2]*va[2]);

    // Computing Lift and Drag   

    // computing t1 

    t1[0] = sin(*theta)*cos(*phi);
    t1[1] = sin(*phi);
    t1[2] = cos(*theta)*cos(*phi);

    // Computing t2 = t1 x vrel

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
    printf("\ni=%d\n", it);
    printf("va[0]=%f, va[1]=%f, va[2]=%f\n", va[0], va[1], va[2]);

    L[0] = *lc*(t2[0]*sin(mu) + t3[0]*cos(mu));
    L[1] = *lc*(t2[1]*sin(mu) + t3[1]*cos(mu));
    L[2] = *lc*(t2[2]*sin(mu) + t3[2]*cos(mu));
    printf("L[0]=%f, L[1]=%f, L[2]=%f\n", L[0], L[1], L[2]);

    *l0 = L[0];
    *l1 = L[1];
    *l2 = L[2];

    L_mod = sqrt(L[0]*L[0] + L[1]*L[1] + L[2]*L[2]); 

    *dc = 0.5*rho*CD_alpha[alpha]*A*Va_mod;

    D[0] = -*dc*va[0];
    D[1] = -*dc*va[1];
    D[2] = -*dc*va[2];
    printf("D[0]=%f, D[1]=%f, D[2]=%f\n", D[0], D[1], D[2]);

    *d0 = D[0];
    *d1 = D[1];
    *d2 = D[2];

    *dc = 0.5*rho*CD_alpha[alpha]*A*Va_mod*sqrt(va[0]*va[0] + va[1]*va[1] + va[2]*va[2]);

    F_aer[0] = L[0] + D[0];
    F_aer[1] = L[1] + D[1];
    F_aer[2] = L[2] + D[2];



    // HOUSKA

    // EFFECTIVE WIND IN THE KITE`S SYSTEM :

    Psi = 0;

    we[0] = -va[0];
    we[1] = -va[1];
    we[2] = -va[2];

    wep[0] = we[0] - (t1[0]*we[0] + t1[1]*we[1] + t1[2]*we[2])*t1[0];
    wep[1] = we[1] - (t1[0]*we[0] + t1[1]*we[1] + t1[2]*we[2])*t1[1];
    wep[2] = we[2] - (t1[0]*we[0] + t1[1]*we[1] + t1[2]*we[2])*t1[2];

    wep_mod = sqrt(wep[0]*wep[0] + wep[1]*wep[1] + wep[2]*wep[2]);

    ew[0] = wep[0]/wep_mod;
    ew[1] = wep[1]/wep_mod;
    ew[2] = wep[2]/wep_mod;

    e0[0] = t1[1]*ew[2] - t1[2]*ew[1];
    e0[1] = t1[2]*ew[0] - t1[0]*ew[2];
    e0[2] = t1[0]*ew[1] - t1[1]*ew[0];

    scal_norm = (t1[0]*we[0] + t1[1]*we[1] + t1[2]*we[2])/wep_mod;
    mod_sqrt = (t1[0]*we[0] + t1[1]*we[1] + t1[2]*we[2])*(t1[0]*we[0] + t1[1]*we[1] + t1[2]*we[2])/(wep_mod*wep_mod);

    et[0] = sin(Psi)*t1[0] - scal_norm*sin(Psi)*ew[0] + 
        sqrt( cos(Psi)*cos(Psi) - mod_sqrt*sin(Psi)*sin(Psi)*e0[0]  );

    et[1] = sin(Psi)*t1[1] - scal_norm*sin(Psi)*ew[1] + 
        sqrt( cos(Psi)*cos(Psi) - mod_sqrt*sin(Psi)*sin(Psi)*e0[1]  );

    et[2] = sin(Psi)*t1[2] - scal_norm*sin(Psi)*ew[2] + 
        sqrt( cos(Psi)*cos(Psi) - mod_sqrt*sin(Psi)*sin(Psi)*e0[2]  );

    printf("et[0]=%f, et[1]=%f, et[2]=%f\n", et[0], et[1], et[2]);

    // CALCULATION OF THE KITE`S TRANSVERSAL AXIS :
    // ---------------------------------------------------------------

    //nwep =  sqrt(we[1]*we[1] + we[2]*we[2]);
    //nwe  =  sqrt(we[0]*we[0] + we[1]*we[1] + we[2]*we[2]);
    //eta  =  asin( we[0]*tan(Psi)/ nwep );
    //printf("newp=%f, new=%f, eta=%f\n", nwep, nwe, eta);

    // ---------------------------------------------------------------

    //ewep[0] = 0.0;
    //ewep[1] = we[1] / nwep;
    //ewep[2] = we[2] / nwep;

    // ---------------------------------------------------------------

    //et[0] = sin(Psi);
    //et[1] = (-cos(Psi)*sin(eta))*ewep[1] - (cos(Psi)*cos(eta))*ewep[2];
    //et[2] = (-cos(Psi)*sin(eta))*ewep[2] + (cos(Psi)*cos(eta))*ewep[1];
    //printf("et[0] = %f, et[1]=%f, et[2]=%f\n", et[0], et[1], et[2]);

    // SUM OF GRAVITATIONAL AND LIFTING FORCE :
    // ---------------------------------------------------------------

    /*Fg[0] =  m*g*cos(theta);
    Fg[1] =  m*g*0.0;
    Fg[2] =  m*g*sin(theta);*/


    // SUM OF THE AERODYNAMIC FORCES :
    // ---------------------------------------------------------------

    Caer = 0.5*rho*A*Va_mod;
    F_aerh[0] =  Caer*( CL_alpha[alpha]*(we[1]*et[2]-we[2]*et[1]) + CD_alpha[alpha]*we[0]);
    F_aerh[1] =  Caer*( CL_alpha[alpha]*(we[2]*et[0]-we[0]*et[2]) + CD_alpha[alpha]*we[1]);
    F_aerh[2] =  Caer*( CL_alpha[alpha]*(we[0]*et[1]-we[1]*et[0]) + CD_alpha[alpha]*we[2]);

    L1[0] = Caer*CL_alpha[alpha]*(we[1]*et[2]-we[2]*et[1]);
    L1[1] = Caer*CL_alpha[alpha]*(we[2]*et[0]-we[0]*et[2]);
    L1[2] = Caer*CL_alpha[alpha]*(we[0]*et[1]-we[1]*et[0]);
    printf("L1[0]=%f, L1[1]=%f, L1[2]=%f\n", L1[0], L1[1], L1[2]);

    D1[0] = Caer*CD_alpha[alpha]*we[0];
    D1[1] = Caer*CD_alpha[alpha]*we[1];
    D1[2] = Caer*CD_alpha[alpha]*we[2];
    printf("D1[0]=%f, D1[1]=%f, D1[2]=%f\n", D1[0], D1[1], D1[2]);

    // Compute tension

    // ===================== CASE 1) BLOCK NOT MOVING ( |v-block| < 10E-6 ) ==================

    v_block_mod = sqrt(v_block[0]*v_block[0] + v_block[1]*v_block[1]);

    if ( v_block_mod < V_THRESHOLD ){

        Txy_norm = sqrt(cos(*phi)*cos(*phi)*sin(*theta)*sin(*theta) + sin(*phi)*sin(*phi));

        // |Mg| > |Tz|

        denom1 = R*(m+m_block)/(m*m_block)
            - cos(*theta)*cos(*phi)/m_block*( r_diff[2] - coeff_friction*(cos(*phi)*sin(*theta)*r_diff[0] + sin(*phi)*r_diff[1])/Txy_norm );

        T1 = (F_aer[0]*r_diff[0] + F_aer[1]*r_diff[1] + F_aer[2]*r_diff[2])/m
        + (v_diff[0]*v_diff[0] + v_diff[1]*v_diff[1] + v_diff[2]*v_diff[2])
        - g*(r_diff[2] - coeff_friction*(cos(*phi)*sin(*theta)*r_diff[0] + sin(*phi)*r_diff[1])/Txy_norm );

        T1 = T1/denom1;

        // |Mg| < |Tz|

        denom2 = R*(m+m_block)/(m*m_block)
            - cos(*theta)*cos(*phi)/m_block*( r_diff[2] + coeff_friction*(cos(*phi)*sin(*theta)*r_diff[0] + sin(*phi)*r_diff[1])/Txy_norm );

        T2 = (F_aer[0]*r_diff[0] + F_aer[1]*r_diff[1] + F_aer[2]*r_diff[2])/m
        + (v_diff[0]*v_diff[0] + v_diff[1]*v_diff[1] + v_diff[2]*v_diff[2])
        - g*(r_diff[2] + coeff_friction*(cos(*phi)*sin(*theta)*r_diff[0] + sin(*phi)*r_diff[1])/Txy_norm );

        T2 = T2/denom2;

        if ( m_block*g > T1*cos(*theta)*cos(*phi) ){
            *sector = 1;
            *T = T1;
        } else if ( m_block*g <= T2*cos(*theta)*cos(*phi) ){
            *sector = 2;
            *T = T2;
        } else { 
            printf("ERROR!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"); 
        }

        Tension[0] = *T*cos(*phi)*sin(*theta);
        Tension[1] = *T*sin(*phi);
        Tension[2] = *T*cos(*theta)*cos(*phi);

        N = m_block*g - Tension[2];

        F_friction[0] = -coeff_friction*fabs(N)*cos(*phi)*sin(*theta)/Txy_norm;
        F_friction[1] = -coeff_friction*fabs(N)*sin(*phi)/Txy_norm;

        *F_attr = sqrt(F_friction[0]*F_friction[0] + F_friction[1]*F_friction[1]);

        a_block[0] = (Tension[0] + F_friction[0])/m_block;
        a_block[1] = (Tension[1] + F_friction[1])/m_block;

        // ============> If the computed tension is bigger than friction force, block moves (so we have to do nothing).

        // ============> If not, recompute friction as: F_friction = -Tension_xy; which gives a_block[0] = 0

        if ( sqrt(Tension[0]*Tension[0] + Tension[1]*Tension[1]) <  *F_attr ){

            *sector = 3;
            
            denom = R*(m+m_block)/(m*m_block) 
            -cos(*phi)*sin(*theta)*r_diff[0]/m_block - sin(*phi)*r_diff[1]/m_block - cos(*phi)*cos(*theta)*r_diff[2]/m_block;
            //- sin(*theta)/m_block*(cos(*phi)*r_diff[0] + sin(*phi)*r_diff[1]) - cos(*theta)/m_block*r_diff[2];

            *T = (F_aer[0]*r_diff[0] + F_aer[1]*r_diff[1] + F_aer[2]*r_diff[2])/m
            + (v_diff[0]*v_diff[0] + v_diff[1]*v_diff[1] + v_diff[2]*v_diff[2]) - g*r_diff[2];

            *T = *T/denom;

            Tension[0] = *T*cos(*phi)*sin(*theta);
            Tension[1] = *T*sin(*phi);
            Tension[2] = *T*cos(*theta)*cos(*phi);

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

        denom1 = R*(m+m_block)/(m*m_block) - cos(*theta)*cos(*phi)/m_block*(r_diff[2] - 
            coeff_friction/v_block_mod*(v_block[0]*r_diff[0] + v_block[1]*r_diff[1]) );

        T1 = (F_aer[0]*r_diff[0] + F_aer[1]*r_diff[1] + F_aer[2]*r_diff[2])/m
            + (v_diff[0]*v_diff[0] + v_diff[1]*v_diff[1] + v_diff[2]*v_diff[2])
            - g*(r_diff[2] - coeff_friction/v_block_mod*(v_block[0]*r_diff[0] + v_block[1]*r_diff[1]) );

        T1 = T1/denom1;

        // |Mg| < |Tz|

        denom2 = R*(m+m_block)/(m*m_block) - cos(*theta)*cos(*phi)/m_block*(r_diff[2] + 
            coeff_friction/v_block_mod*(v_block[0]*r_diff[0] + v_block[1]*r_diff[1]) );

        T2 = (F_aer[0]*r_diff[0] + F_aer[1]*r_diff[1] + F_aer[2]*r_diff[2])/m
            + (v_diff[0]*v_diff[0] + v_diff[1]*v_diff[1] + v_diff[2]*v_diff[2])
            - g*(r_diff[2] + coeff_friction/v_block_mod*(v_block[0]*r_diff[0] + v_block[1]*r_diff[1]) );

        T2 = T2/denom2;

        if ( m_block*g > T1*cos(*theta)*cos(*phi) ) {
            *sector = 4;
            *T = T1;
        } else if ( m_block*g <= T2*cos(*theta)*cos(*phi) ){
            *sector = 5;
            *T = T2;
        } else {
            printf("ERROR!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"); 
        }

        Tension[0] = *T*cos(*phi)*sin(*theta);
        Tension[1] = *T*sin(*phi);
        Tension[2] = *T*cos(*theta)*cos(*phi);

        N = m_block*g - Tension[2];

        F_friction[0] = -coeff_friction*fabs(N)*v_block[0]/v_block_mod;
        F_friction[1] = -coeff_friction*fabs(N)*v_block[1]/v_block_mod;

        *F_attr = sqrt(F_friction[0]*F_friction[0] + F_friction[1]*F_friction[1]);

        a_block[0] = ( Tension[0] + F_friction[0] )/m_block;
        a_block[1] = ( Tension[1] + F_friction[1] )/m_block;

    }

    a_block[2] = 0;

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

}

#endif

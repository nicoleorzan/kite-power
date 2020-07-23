#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "constants1.h"

#ifndef __dynamics__
#define __dynamics__

#define CONST 10000

//#define DEBUG

double mod_we;
double we[3];
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

double er[3];
double v_block_mod;

double wep_mod;
double scal_norm, mod_sqrt;
double wep[3];
double ew[3], e0[3];
double et[3];
double Caer;
double sqrt_var;

double nwep;
double ewep[3];

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
                            int alpha, double Psi,
                            double * W, double * lc, double * dc,
                            double * T, double *F_attr, int it, int * sector, 
                            double *l0, double *l1, double *l2, double *d0, double *d1, double *d2){

    //printf("\ni=%d, houska\n", it);                              
    
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

    we[0] = W[0] - vk[0];              // Effective wind on x
    we[1] = W[1] - vk[1];              // Effective wind on y
    we[2] = W[2] - vk[2];              // Effective wind on z
    //printf("we[0]=%f, we[1]=%f, we[2]=%f\n", we[0], we[1], we[2]);

    mod_we = sqrt(we[0]*we[0] + we[1]*we[1] + we[2]*we[2]);
    //printf("we[0]=%f, we[1]=%f, we[2]=%f\n", we[0]/mod_we, we[1]/mod_we, we[2]/mod_we);

    // Computing Lift and Drag   

    // computing er

    er[0] = sin(*theta)*cos(*phi);
    er[1] = sin(*theta)*sin(*phi);
    er[2] = cos(*theta);
    //printf("er[0]=%f, er[1]=%f, er[2]=%f\n", er[0], er[1], er[2]);

    // calcolo la componente del vento effettivo sul piano theta/phi sottraendo ad ogni componente il prod
    // scalare tra We ed er:
    // We_piano = We - (er*we)*er

    wep[0] = we[0] - (er[0]*we[0] + er[1]*we[1] + er[2]*we[2])*er[0];
    wep[1] = we[1] - (er[0]*we[0] + er[1]*we[1] + er[2]*we[2])*er[1];
    wep[2] = we[2] - (er[0]*we[0] + er[1]*we[1] + er[2]*we[2])*er[2];
    //printf("wep[0]=%f, wep[1]=%f, wep[2]=%f\n", wep[0], wep[1], wep[2]);

    // normalizzo wep ==> ew

    wep_mod = sqrt(wep[0]*wep[0] + wep[1]*wep[1] + wep[2]*wep[2]);

    ew[0] = wep[0]/wep_mod;
    ew[1] = wep[1]/wep_mod;
    ew[2] = wep[2]/wep_mod;
    //printf("ew[0]=%f, ew[1]=%f, ew[2]=%f\n", ew[0], ew[1], ew[2]);

    // calcolo e0, vettore perpendicolare ad er ed ew

    e0[0] = er[1]*ew[2] - er[2]*ew[1];
    e0[1] = er[2]*ew[0] - er[0]*ew[2];
    e0[2] = er[0]*ew[1] - er[1]*ew[0];
    //printf("e0[0]=%f, e0[1]=%f, e0[2]=%f\n", e0[0], e0[1], e0[2]);

    // per calcolare et, lo scompongo sulla terna er, ew, e0
    // et = ar*er + aw*ew + a0*e0 ==> vedi eq 8 del paper che ti ho mandato

    scal_norm = (er[0]*we[0] + er[1]*we[1] + er[2]*we[2])/wep_mod;
    mod_sqrt = (er[0]*we[0] + er[1]*we[1] + er[2]*we[2])*(er[0]*we[0] + er[1]*we[1] + er[2]*we[2])/(wep_mod*wep_mod);
    //printf("(er[0]*we[0] + er[1]*we[1] + er[2]*we[2])=%f\n",  (er[0]*we[0] + er[1]*we[1] + er[2]*we[2]));

    sqrt_var = cos(Psi)*cos(Psi) - mod_sqrt*sin(Psi)*sin(Psi);

    /*if ( sqrt_var < 0){
        sqrt_var = -sqrt_var;
    }*/
    //printf("sqrt(sqrt_var)=%f\n", sqrt(sqrt_var));

    et[0] = sin(Psi)*er[0] - scal_norm*sin(Psi)*ew[0] + 
        sqrt( sqrt_var )*e0[0];
        //sqrt( cos(Psi)*cos(Psi) - mod_sqrt*sin(Psi)*sin(Psi) )*e0[0];

    et[1] = sin(Psi)*er[1] - scal_norm*sin(Psi)*ew[1] + 
        sqrt( sqrt_var )*e0[1];
        //sqrt( cos(Psi)*cos(Psi) - mod_sqrt*sin(Psi)*sin(Psi) )*e0[1];

    et[2] = sin(Psi)*er[2] - scal_norm*sin(Psi)*ew[2] + 
        sqrt( sqrt_var )*e0[2];
    printf("et[0]=%f, et[1]=%f, et[2]=%f\n", et[0], et[1], et[2]);

    Caer = 0.5*rho*A*mod_we;

    // scrivo lift e drag come L = C*(we x et), D = C*el. non uso we normalizzato perche` 
    // cosi` ho dentro anche la componente della velocita che mi serve per calcoalre il lift

    L[0] = Caer*CL_alpha[alpha]*(we[1]*et[2]-we[2]*et[1]);
    L[1] = Caer*CL_alpha[alpha]*(we[2]*et[0]-we[0]*et[2]);
    L[2] = Caer*CL_alpha[alpha]*(we[0]*et[1]-we[1]*et[0]);
    //printf("L[0]=%f, L[1]=%f, L[2]=%f\n", L[0], L[1], L[2]);
    
    *lc =  Caer*CL_alpha[alpha]*mod_we;

    *l0 = L[0];
    *l1 = L[1];
    *l2 = L[2];

    D[0] = Caer*CD_alpha[alpha]*we[0];
    D[1] = Caer*CD_alpha[alpha]*we[1];
    D[2] = Caer*CD_alpha[alpha]*we[2];
    //printf("D[0]=%f, D[1]=%f, D[2]=%f\n\n", D[0], D[1], D[2]);

    *d0 = D[0];
    *d1 = D[1];
    *d2 = D[2];

    *dc = Caer*CD_alpha[alpha]*mod_we;

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

    // ========================== CASE 2) BLOCK MOVING (|v| >= 10E-6) ===> Fmu = -Psi*|N|*vx/|vx| =====================

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

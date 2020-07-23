#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "constants1.h"

#ifndef __dynamics__
#define __dynamics__

//#define DEBUG

double Va_mod;
double va[2];
double L[2];
double D[2];
double beta;
double Fg[2] = {0, -m*g};
double F_friction;
double N;
double detA, detA0, detA1, detA2;
double csi;
double C[3];
double Tension[2];
double Ftot[2]; 
double T1, a_block1, ddtheta1;
double T2, a_block2, ddtheta2;
double t2;
double t3[2];

//double theta_star;

void variables_initialization(double * rk, double * vk, double * ak, 
                            double _theta, double _vtheta,
                            double * r_block, double * v_block, double * a_block,
                            double * theta){ // Angle, velocity and acceleration values
    r_block[0] = 0;
    r_block[1] = 0;

    v_block[0] = 0;
    v_block[1] = 0;

    a_block[0] = 0;
    a_block[1] = 0;

    theta[0] = _theta;
    theta[1] = _vtheta;
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
                            int alpha,
                            double * W, double * lc, double * dc,
                            double * T, double *F_attr, int it, int * sector, int* et_val){

    va[0] = vk[0] - W[0];              // Apparent velocity on x
    va[1] = vk[1] - W[1];              // Apparent velocity on z

    Va_mod = sqrt(va[0]*va[0] + va[1]*va[1]);

    // Computing Lift and Drag

    beta = atan2(va[1], va[0]);

    *lc = 0.5*rho*CL_alpha[alpha]*A*Va_mod*Va_mod;

    *dc = 0.5*rho*CD_alpha[alpha]*A*Va_mod*Va_mod;

    //t2 = (cos(*theta)*va[1] - sin(*theta)*va[0])/fabs(cos(*theta)*va[1] - sin(*theta)*va[0]);
    if (it == 0){
        t2 = (cos(*theta)*va[1] - sin(*theta)*va[0])/fabs(cos(*theta)*va[1] - sin(*theta)*va[0]);
        *et_val = t2;
    }

    t3[0] = va[1]*t2/Va_mod;
    t3[1] = -va[0]*t2/Va_mod;

    L[0] = *lc*t3[0];
    L[1] = *lc*t3[1];

    D[0] = *dc*cos(beta + PI);
    D[1] = *dc*sin(beta + PI);

    // ===================== CASE 1) BLOCK NOT MOVING ( |v-block| < 10E-6 ) ==================

    if ( fabs(v_block[0]) < V_THRESHOLD ){

        // |Mg| > |Tz|

        csi = - cos(theta[0]) - coeff_friction*sin(theta[0])*cos(theta[0])/fabs(cos(*theta));

        detA = m*(m*R*cos(theta[0])*csi) 
            - m*R*sin(theta[0])*m_block*sin(theta[0]) 
            - cos(theta[0])*m*m_block*R*cos(theta[0]);

        C[0] = D[0] + L[0] + m*R*theta[1]*theta[1]*cos(theta[0]);
        C[1] = D[1] + L[1] - m*g + m*R*theta[1]*theta[1]*sin(theta[0]);
        C[2] = -coeff_friction*m_block*g*cos(theta[0])/fabs(cos(*theta));

        detA0 = C[0]*m*R*cos(theta[0])*csi 
            + m*R*sin(theta[0])*(C[1]*csi - C[2]*sin(theta[0]))
            - cos(theta[0])*C[2]*m*R*cos(theta[0]);

        detA1 = m*(C[1]*csi - C[2]*sin(theta[0]))
            + C[0]*m_block*sin(theta[0])
            - cos(theta[0])*m_block*C[1];

        detA2 = m*m*R*cos(theta[0])*C[2]
            - m*R*sin(theta[0])*m_block*C[1]
            - C[0]*m*m_block*R*cos(theta[0]);

        a_block1 = detA0/detA;
        ddtheta1 = detA1/detA;
        T1 = detA2/detA;

        // |Mg| < |Tz|

        csi = - cos(theta[0]) + coeff_friction*sin(theta[0])*cos(theta[0])/fabs(cos(*theta));

        detA = m*(m*R*cos(theta[0])*csi) 
            - m*R*sin(theta[0])*m_block*sin(theta[0]) 
            - cos(theta[0])*m*m_block*R*cos(theta[0]);

        C[0] = D[0] + L[0] + m*R*theta[1]*theta[1]*cos(theta[0]);
        C[1] = D[1] + L[1] - m*g + m*R*theta[1]*theta[1]*sin(theta[0]);
        C[2] = coeff_friction*m_block*g*cos(theta[0])/fabs(cos(*theta));

        detA0 = C[0]*m*R*cos(theta[0])*csi 
            + m*R*sin(theta[0])*(C[1]*csi - C[2]*sin(theta[0]))
            - cos(theta[0])*C[2]*m*R*cos(theta[0]);

        detA1 = m*(C[1]*csi - C[2]*sin(theta[0]))
            + C[0]*m_block*sin(theta[0])
            - cos(theta[0])*m_block*C[1];

        detA2 = m*m*R*cos(theta[0])*C[2]
            - m*R*sin(theta[0])*m_block*C[1]
            - C[0]*m*m_block*R*cos(theta[0]);

        a_block2 = detA0/detA;
        ddtheta2 = detA1/detA;
        T2 = detA2/detA;

        if ( m_block*g > T1*sin(*theta) ){
            *sector = 1;
            a_block[0] = a_block1;
            theta[2] = ddtheta1;
            *T = T1;
        } else if ( m_block*g < T2*sin(*theta) ){
            *sector = 2;
            a_block[0] = a_block2;
            theta[2] = ddtheta2;
            *T = T2;            
        } else { 
            printf("ERROR DYANAMICS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"); 
        }

        Tension[0] = *T*cos(theta[0]);
        Tension[1] = *T*sin(theta[0]);   

        N = m_block*g - Tension[1];

        F_friction = -coeff_friction*fabs(N)*cos(theta[0]);
        *F_attr = F_friction;

        // ============> If the computed tension is bigger than friction force, block moves (so we have to do nothing).

        // ============> If not, recompute friction as: F_friction = -Tension[0]; which gives a_block[0] = 0

        if ( fabs(Tension[0]) < fabs(F_friction) ){

            *sector = 3;
        
            csi = 0;

            detA = m*(m*R*cos(theta[0])*csi) 
                - m*R*sin(theta[0])*m_block*sin(theta[0]) 
                - cos(theta[0])*m*m_block*R*cos(theta[0]);

            C[0] = D[0] + L[0] + m*R*theta[1]*theta[1]*cos(theta[0]);
            C[1] = D[1] + L[1] - m*g + m*R*theta[1]*theta[1]*sin(theta[0]);
            C[2] = 0;

            detA0 = C[0]*m*R*cos(theta[0])*csi 
                + m*R*sin(theta[0])*(C[1]*csi - C[2]*sin(theta[0]))
                - cos(theta[0])*C[2]*m*R*cos(theta[0]);

            detA1 = m*(C[1]*csi - C[2]*sin(theta[0]))
                + C[0]*m_block*sin(theta[0])
                - cos(theta[0])*m_block*C[1];

            detA2 = m*m*R*cos(theta[0])*C[2]
                - m*R*sin(theta[0])*m_block*C[1]
                - C[0]*m*m_block*R*cos(theta[0]);

            a_block[0] = detA0/detA; // should be zero
            theta[2] = detA1/detA;
            *T = detA2/detA;

            Tension[0] = *T*cos(theta[0]);
            Tension[1] = *T*sin(theta[0]);

            N = m_block*g - Tension[1];

            F_friction = -Tension[0];
            *F_attr = F_friction;
        }

    // ========================== CASE 2) BLOCK MOVING (|v| >= 10E-6) ===> Fmu = -mu*|N|*vx/|vx| =====================

    } else {

        // |Mg| > |Tz|

        csi = - cos(theta[0]) - coeff_friction*sin(theta[0])*v_block[0]/fabs(v_block[0]);

        detA = m*(m*R*cos(theta[0])*csi) 
            - m*R*sin(theta[0])*m_block*sin(theta[0]) 
            - cos(theta[0])*m*m_block*R*cos(theta[0]);

        C[0] = D[0] + L[0] + m*R*theta[1]*theta[1]*cos(theta[0]);
        C[1] = D[1] + L[1] - m*g + m*R*theta[1]*theta[1]*sin(theta[0]);
        C[2] = -coeff_friction*m_block*g*v_block[0]/fabs(v_block[0]);

        detA0 = C[0]*m*R*cos(theta[0])*csi 
            + m*R*sin(theta[0])*(C[1]*csi - C[2]*sin(theta[0]))
            - cos(theta[0])*C[2]*m*R*cos(theta[0]);

        detA1 = m*(C[1]*csi - C[2]*sin(theta[0]))
            + C[0]*m_block*sin(theta[0])
            - cos(theta[0])*m_block*C[1];

        detA2 = m*m*R*cos(theta[0])*C[2]
            - m*R*sin(theta[0])*m_block*C[1]
            - C[0]*m*m_block*R*cos(theta[0]);

        a_block1 = detA0/detA;
        ddtheta1 = detA1/detA;
        T1 = detA2/detA;

        //  Mg < Tz 

        csi = - cos(theta[0]) + coeff_friction*sin(theta[0])*v_block[0]/fabs(v_block[0]);

        detA = m*(m*R*cos(theta[0])*csi) 
            - m*R*sin(theta[0])*m_block*sin(theta[0]) 
            - cos(theta[0])*m*m_block*R*cos(theta[0]);

        C[0] = D[0] + L[0] + m*R*theta[1]*theta[1]*cos(theta[0]);
        C[1] = D[1] + L[1] - m*g + m*R*theta[1]*theta[1]*sin(theta[0]);
        C[2] = coeff_friction*m_block*g*v_block[0]/fabs(v_block[0]);

        detA0 = C[0]*m*R*cos(theta[0])*csi 
            + m*R*sin(theta[0])*(C[1]*csi - C[2]*sin(theta[0]))
            - cos(theta[0])*C[2]*m*R*cos(theta[0]);

        detA1 = m*(C[1]*csi - C[2]*sin(theta[0]))
            + C[0]*m_block*sin(theta[0])
            - cos(theta[0])*m_block*C[1];

        detA2 = m*m*R*cos(theta[0])*C[2]
            - m*R*sin(theta[0])*m_block*C[1]
            - C[0]*m*m_block*R*cos(theta[0]);

        a_block2 = detA0/detA;
        ddtheta2 = detA1/detA;
        T2 = detA2/detA;

        if ( m_block*g > T1*sin(*theta) ){
            *sector = 4;
            a_block[0] = a_block1;
            theta[2] = ddtheta1;
            *T = T1;
        } else if ( m_block*g < T2*sin(*theta) ){
            *sector = 5;
            a_block[0] = a_block2;
            theta[2] = ddtheta2;
            *T = T2;            
        } else { 
            printf("ERROR DYANAMICS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"); 
        }

        Tension[0] = *T*cos(theta[0]);
        Tension[1] = *T*sin(theta[0]);    

        N = m_block*g - Tension[1];

        F_friction = -coeff_friction*fabs(N)*v_block[0]/fabs(v_block[0]);
        *F_attr = F_friction;
    }              

    //a_block[1] = 0;

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

    //Ftot[0] = L[0] + D[0] + Fg[0] - Tension[0];
    //Ftot[1] = L[1] + D[1] + Fg[1] - Tension[1];

}

#endif

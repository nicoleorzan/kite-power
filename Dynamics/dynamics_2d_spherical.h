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
#define coeff_friction 0.4      // Block friction coefficient
#define m_blocco 100.0          // Block mass

double W;
double Vw[2], Va_mod;
double Va[2];
double L[2], Lc;
double D[2], Dc;
double Fg[2];
//double gamma_;                // gliding angle
double F_perpendic;
double theta_star;
double beta;
double T_orizz;
double Tension[2];
double Ftot[2]; 

void variables_initialization(double * x_blocco, double * v_blocco, double * a_blocco,
                              double * theta,              // Angle, velocity and acceleration values
                              double * rk, double * vk, double * ak, 
                              double theta0, double vtheta0){
    x_blocco[0] = 0;
    x_blocco[1] = 0;
    v_blocco[0] = 0;
    v_blocco[1] = 0;
    a_blocco[0] = 0;
    a_blocco[1] = 0;

    theta[0] = theta0;
    theta[1] = vtheta0;
    //theta[0] = ((float)rand()/(float)(RAND_MAX)) * 3.14;    // Random num between 0 and 3.14
    //theta[1] = ((float)rand()/(float)(RAND_MAX)) - 0.5;     // Random num between -0.5 and 0.5
    theta[2] = 0.;

    rk[0] = R*cos(theta[0]) + x_blocco[0];
    rk[1] = R*sin(theta[0]);
    vk[0] = v_blocco[0] - R*sin(theta[0])*theta[1];
    vk[1] = R*cos(theta[0])*theta[1];
    ak[0] = a_blocco[0] - R*cos(theta[0])*theta[1]*theta[1] - R*sin(theta[0])*theta[2];
    ak[1] = -R*sin(theta[0])*theta[1]*theta[1] + R*cos(theta[0])*theta[2];
}

void integration_trajectory(double * rk, double * vk, double * ak, // Kite variables
                            double * x_blocco, double * v_blocco, double * a_blocco, // Block variables
                            double * theta, // Angle, velocity and acceleration values
                            double * T, double * power, double Cl, double Cd,
                            double Wx, double Wy, double * lift, double * drag){

Vw[0] = Wx;                          // Wind velocity on x
Vw[1] = Wy;                          // Wind velocity on z

Va[0] = vk[0] - Vw[0];              // Apparent velocity on x
Va[1] = vk[1] - Vw[1];              // Apparent velocity on z

Va_mod = sqrt(Va[0]*Va[0] + Va[1]*Va[1]);

//gamma_ = atan(Va[1]/Va[0]);

beta = atan2(Va[1], Va[0]);

Lc = 0.5*rho*Cl*A*Va_mod*Va_mod;
*lift = Lc;

Dc = 0.5*rho*Cd*A*Va_mod*Va_mod;
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

theta_star = atan((Lc - m*g)/Dc);

// Gravitional force Kite (x, z)

Fg[0] = 0;                      
Fg[1] = -m*g;                   

// Theta acceleration

theta[2] = (a_blocco[0]*sin(theta[0]*m) - m*g*cos(theta[0]) - L[0]*sin(theta[0]) - D[0]*sin(theta[0]) + L[1]*cos(theta[0]) + D[1]*cos(theta[0]) )/(m*R);

// Tension next step

*T = - a_blocco[0]*cos(theta[0])*m + m*R*theta[1]*theta[1] - m*g*sin(theta[0]) + L[0]*cos(theta[0]) + L[1]*sin(theta[0]) + D[0]*cos(theta[0]) + D[1]*sin(theta[0]);

// Compute Block motion

F_perpendic = coeff_friction*fabs(m_blocco*g - *T*sin(theta[0]));

if ( abs(v_blocco[0]) < 10E-6 ){
    if ( abs(*T*cos(theta[0])) > abs(F_perpendic)  ){
        if (cos(theta[0]) > 0){
            a_blocco[0] = (*T*cos(theta[0]) - F_perpendic )/m_blocco;
        }
        else{
            a_blocco[0] = (*T*cos(theta[0]) + F_perpendic )/m_blocco;
        }
    }
    else { a_blocco[0] = 0; } 
}
else if ( v_blocco[0] > 0 ){
    a_blocco[0] = (*T*cos(theta[0]) - F_perpendic )/m_blocco;
}
else if ( v_blocco[0] < 0 ){
    a_blocco[0] = (*T*cos(theta[0]) + F_perpendic )/m_blocco;
}

v_blocco[0] = v_blocco[0] + h*a_blocco[0];

// Case 1: First quadrant, 0 < theta < PI/2

/*if (cos(theta[0]) > 0){

    if (*T*cos(theta[0]) > F_perpendic ){
        a_blocco[0] = (*T*cos(theta[0]) - F_perpendic )/m_blocco; // deve essere positiva
        v_blocco[0] = v_blocco[0] + h*a_blocco[0];
    }
    else if (*T*cos(theta[0]) < F_perpendic && v_blocco[0] > 0.){
        a_blocco[0] = (*T*cos(theta[0]) - F_perpendic )/m_blocco; // deve essere negativa
        v_blocco[0] = v_blocco[0] + h*a_blocco[0];
    }
    else if ( *T*cos(theta[0]) <= F_perpendic && v_blocco[0] <= 0.){ //se blocco fermo e tensione piu debole della f attrito
        a_blocco[0] = 0;
        v_blocco[0] = 0.;
    }

} 

// Case 2: Second quadrant, PI/2 < theta < PI

else if (cos(theta[0]) < 0){

    if (-*T*cos(theta[0]) > F_perpendic ){
        a_blocco[0] = (*T*cos(theta[0]) + F_perpendic )/m_blocco; // deve essere negativa
        v_blocco[0] = v_blocco[0] + h*a_blocco[0];
    }
    else if (-*T*cos(theta[0]) < F_perpendic && v_blocco[0] < 0.){
        a_blocco[0] = (*T*cos(theta[0]) + F_perpendic )/m_blocco; // deve essere positiva
        v_blocco[0] = v_blocco[0] + h*a_blocco[0];
    }
    else if ( -*T*cos(theta[0]) <= F_perpendic && v_blocco[0] >= 0.){ //se blocco fermo e tensione piu debole della f attrito
        a_blocco[0] = 0;
        v_blocco[0] = 0.;
    }

}*/

x_blocco[0] = x_blocco[0] + h*v_blocco[0]; 

theta[1] = theta[1] + h*theta[2];

theta[0] = theta[0] + h*theta[1];

// Kite position, velocity and acceleration update (x, z)

rk[0] = x_blocco[0] + R*cos(theta[0]);

rk[1] = R*sin(theta[0]);

vk[0] = v_blocco[0] - R*sin(theta[0])*theta[1];

vk[1] = R*cos(theta[0])*theta[1];

ak[0] = a_blocco[0] - R*cos(theta[0])*theta[1]*theta[1] - R*sin(theta[0])*theta[2];

ak[1] = -R*sin(theta[0])*theta[1]*theta[1] + R*cos(theta[0])*theta[2];

// Check: Sum of total forces

Tension[0] = -*T*cos(theta[0]);

Tension[1] = -*T*sin(theta[0]);

Ftot[0] = L[0] + D[0] + Fg[0] + Tension[0];

Ftot[1] = L[1] + D[1] + Fg[1] + Tension[1];

// Power

//*power = *T*dpos[0]*sin(pos[1]); // ?

#ifdef DEBUG
    printf("thetastar = %f\n", theta_star);
    printf("beta=%f\n", beta);
    printf("Vw[0]= %f, Vw[1]=%f\n", Vw[0], Vw[1]);
    printf("Vkx=%f, Vkz=%f\n", vk[0], vk[1]); 
    printf("Va_mod=%f, Va[0]=%f, Va[1]=%f\n", Va_mod, Va[0], Va[1]); 
    //printf("gamma_=%f\n", gamma_);
    printf("L[0]=%f\n", L[0]); 
    printf("L[1]=%f\n", L[1]); 
    printf("D[0]=%f\n", D[0]); 
    printf("D[1]=%f\n", D[1]); 
    printf("Fg[0]=%f\n", Fg[0]); 
    printf("Fg[1]=%f\n", Fg[1]); 
    printf("T=%f, a_blocco[0]=%f, acc_theta = %f\n", *T, a_blocco[0], theta[2]);
    printf("Tx=%f, Ty=%f\n", *T*cos(theta[0]), *T*sin(theta[0]));
    printf("theta[0] = %f, theta[1] = %f\n", theta[0], theta[1]);
    printf("a_blocco[0] = %f, a_blocco[1] = %f\n", a_blocco[0], a_blocco[1]);
    printf("v_blocco[0] = %f, v_blocco[1] = %f\n", v_blocco[0], v_blocco[1]);
    printf("x_blocco[0] = %f, x_blocco[1] = %f\n", x_blocco[0], x_blocco[1]);
    printf("rk[0]= %f, rk[1]=%f\n", rk[0], rk[1]);
    printf("vk[0]= %f, vk[1]=%f\n", vk[0], vk[1]);
    printf("ak[0]= %f, ak[1]=%f\n", ak[0], ak[1]);
    printf("forze totali x:%f\n", Ftot[0]);
    printf("forze totali x:%f\n", Ftot[1]);
    //printf("power %f\n", *power);
    printf("\n");
#endif    

}

#endif

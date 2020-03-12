#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "constants.h"

#ifndef __dynamics__
#define __dynamics__

#define DEBUG

#define h 0.0001
// fisso alpha = 14, l'aereo guarda verso -x
#define Cl 1.0
#define Cd 0.13
#define W 10
//#define coeff_friction 0.4 
//#define m_blocco 10000000

double Vw[2], Va_mod;
double Va[2];
double L[2], Lc;
double D[2], Dc;
double Fg[2];
double gamma_;
double detA, detA1, detA2, detA3;
double c[3];
double theta_star;
double beta;
double Tension[2];
double forzetot[2]; 
double Lr, Ltheta, Dr, Dtheta;

void integration_trajectory(double * rk, double * vk, double * ak, // Kite
                            double * x_blocco, double * v_blocco, double * acc_blocco, // Blocco
                            double * theta, // theta0=angolo, theta1=vel, theta2=acc
                            double * T, double * power){

// W e` il modulo della velocita del vento
Vw[0] = W; // velocita del vento in x
Vw[1] = 0; // velocita del vento in z

Va[0] = vk[0] - Vw[0]; //velocita apparente lungo x
Va[1] = vk[1] - Vw[1]; //velocita apparente lungo z

Va_mod = sqrt(Va[0]*Va[0] + Va[1]*Va[1]);

// Compute gamma_, inclinazione di Va rispetto all'asse orizzontale (gliding angle)

gamma_ = atan(Va[1]/Va[0]); //atan2(Va[1], Va[0]);
beta = atan2(Va[1], Va[0]);
printf("beta=%f\n", beta);
printf("Vw[0]= %f, Vw[1]=%f\n", Vw[0], Vw[1]);
printf("Vkx=%f, Vkz=%f\n", vk[0], vk[1]); 
printf("Va_mod=%f, Va[0]=%f, Va[1]=%f\n", Va_mod, Va[0], Va[1]); 
printf("gamma_=%f\n", gamma_);

Lc = 0.5*rho*Cl*A*Va_mod*Va_mod;

Dc = 0.5*rho*Cd*A*Va_mod*Va_mod;

// Lift

if (beta > - PI/2. && beta < PI/2.){
    printf("if\n");
    L[0] = Lc*cos(beta - PI/2.);
    L[1] = -Lc*sin(beta - PI/2.);
} else { 
    L[0] = Lc*cos(beta - PI/2.); //fabs(Lc*sin(gamma_));
    L[1] = Lc*sin(beta - PI/2.); //fabs(Lc*cos(gamma_));
}

// Drag

D[0] = Dc*cos(beta + PI); //fabs(Dc*cos(gamma_)); 
D[1] = Dc*sin(beta + PI); //fabs(Dc*sin(gamma_));

printf("beta+pi %f\n", beta + PI);
printf("cos(beta+pi) %f\n", cos(beta + PI));
printf("sin(beta+pi) %f\n", sin(beta + PI));

theta_star = atan((Lc - m*g)/Dc);

printf("thetastar = %f\n", theta_star);

// Gravitional force kite

Fg[0] = 0; // forza grav lungo x
Fg[1] = -m*g; // forza grav lungo z

// Calcolo matrice

//theta[0] = atan2(rk[1], rk[0] - x_blocco[0]);

theta[2] = (-m*g*cos(theta[0]) - L[0]*sin(theta[0]) - D[0]*sin(theta[0]) + L[1]*cos(theta[0]) + D[1]*cos(theta[0]) )/(m*R);

*T = m*R*theta[1]*theta[1] - m*g*sin(theta[0]) + L[0]*cos(theta[0]) + L[1]*sin(theta[0]) + D[0]*cos(theta[0]) + D[1]*sin(theta[0]);

v_blocco[0] = v_blocco[0] + h*acc_blocco[0];

x_blocco[0] = x_blocco[0] + h*v_blocco[0]; 

theta[1] = theta[1] + h*theta[2];

theta[0] = theta[0] + h*theta[1];

rk[0] = x_blocco[0] + R*cos(theta[0]);

rk[1] = R*sin(theta[0]);

vk[0] = v_blocco[0] - R*sin(theta[0])*theta[1];

vk[1] = R*cos(theta[0])*theta[1];

ak[0] = acc_blocco[0] - R*cos(theta[0])*theta[1]*theta[1] - R*sin(theta[0])*theta[2];

ak[1] = -R*sin(theta[0])*theta[1]*theta[1] + R*cos(theta[0])*theta[2];

// sum of total forces

Tension[0] = -*T*cos(theta[0]);

Tension[1] = -*T*sin(theta[0]);

forzetot[0] = L[0] + D[0] + Fg[0] + Tension[0];

forzetot[1] = L[1] + D[1] + Fg[1] + Tension[1];

// Power

//*power = *T*dpos[0]*sin(pos[1]); // ?

#ifdef DEBUG
    //printf("Vw[0]= %f, Vw[1]=%f\n", Vw[0], Vw[1]);
    //printf("Va_mod=%f, Va[0]=%f, Va[1]=%f\n", Va_mod, Va[0], Va[1]); 
    //printf("gamma_=%f\n", gamma_);
    printf("L[0]=%f\n", L[0]); 
    printf("L[1]=%f\n", L[1]); 
    printf("D[0]=%f\n", D[0]); 
    printf("D[1]=%f\n", D[1]); 
    printf("Fg[0]=%f\n", Fg[0]); 
    printf("Fg[1]=%f\n", Fg[1]); 
    printf("T=%f, acc_blocco[0]=%f, acc_theta = %f\n", *T, acc_blocco[0], theta[2]);
    printf("Tx=%f, Ty=%f\n", *T*cos(theta[0]), *T*sin(theta[0]));
    printf("theta[0] = %f, theta[1] = %f\n", theta[0], theta[1]);
    printf("acc_blocco[0] = %f, acc_blocco[1] = %f\n", acc_blocco[0], acc_blocco[1]);
    printf("v_blocco[0] = %f, v_blocco[1] = %f\n", v_blocco[0], v_blocco[1]);
    printf("x_blocco[0] = %f, x_blocco[1] = %f\n", x_blocco[0], x_blocco[1]);
    printf("rk[0]= %f, rk[1]=%f\n", rk[0], rk[1]);
    printf("vk[0]= %f, vk[1]=%f\n", vk[0], vk[1]);
    printf("ak[0]= %f, ak[1]=%f\n", ak[0], ak[1]);
    printf("forze totali x:%f\n", forzetot[0]);
    printf("forze totali x:%f\n", forzetot[1]);
    //printf("power %f\n", *power);
    printf("\n");
#endif    

}



#endif

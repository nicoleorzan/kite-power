#ifndef constants
#define constants

#define PI 3.1415926535897932384626433

//physical constants
#define rho 1.225           	// kg/m^3
#define g 9.81              	// m/s^2

// kite
#define m 1 //50.0              // kg
#define A 5 //100               // m^2
#define R 50.0                   // cable length m (in the paper varies between 500 and 1000)

// block
#define m_block 40.0           // Block mass
#define coeff_friction 0.4      // Block friction coefficient

/* Paper: Fagiano et al., "High Altitude Wind Energy Generation Using Controlled Power Kites", 2010 */
#define n_alphas 15
double CL_alpha[n_alphas] = {-0.15, -0.05, 0.05, 0.2, 0.35, 0.45, 0.55, 0.65, 0.75, 0.82, 0.9, 1.0, 1.08, 1.1, 1.05};
double CD_alpha[n_alphas] = {0.005, 0.005, 0.001, 0.005, 0.01, 0.02, 0.03, 0.05, 0.07, 0.09, 0.1, 0.13, 0.18, 0.18, 0.21};
double alphas[n_alphas] = {-8, -6, -4, -2, 0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20};

// integration constants
#define h 0.0001                // Integration dt
#define STEPS 2000000
#define PRINTSTEP 100
#define THRESHOLD 10E-5

#define V_THRESHOLD 10E-8

// 2d initial conditions
#define theta0 PI/4.
#define vtheta0 0.

#endif

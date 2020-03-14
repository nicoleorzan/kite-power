#ifndef constants
#define constants

#define PI 3.1415926535897932384626433

//physical constants
#define rho 1.225           	// kg/m^3
#define g 9.81              	// m/s^2

// kite
#define m 10 //50.0              // kg
#define A 25 //100               // m^2
#define R 50.0                   // cable length m (in the paper varies between 500 and 1000)

// block
#define m_block 1000.0           // Block mass

/* Paper: Fagiano et al., "High Altitude Wind Energy Generation Using Controlled Power Kites", 2010 */
#define n_alphas 16
double CL_alpha[n_alphas] = {-0.15, -0.05, 0.05, 0.2, 0.35, 0.42, 0.56, 0.65, 0.77, 0.82, 0.9, 1., 1.3, 1.4, 1.3, 1.1};
double CD_alpha[n_alphas] = {0.005, 0.005, 0.001, 0.005, 0.01, 0.02, 0.03, 0.05, 0.07, 0.09, 0.1, 0.13, 0.18, 0.195, 0.18, 0.21};
double alphas[n_alphas] = {-8, -6, -4, -2, 0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 19, 20};

// integration constants
#define h 0.00001                // Integration dt
#define STEPS 1000000
#define PRINTSTEP 5000

#define V_THRESHOLD 10E-8

#endif

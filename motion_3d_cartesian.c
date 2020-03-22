#include "Dynamics/dynamics_3d_cartesian.h"
//#include "Dynamics/winds.h"
#include <time.h>
#include <string.h>
#include <stdbool.h>

#define _theta0 PI/2.
#define _dtheta0 -1.0
#define phi0 0.
#define dphi0 0.
#define dim 3
#define mu 0.

// ============== FILE INPUT: ATTACK ANGLE AND WIND COEFF AND IF TRAJECTORY IS NEEDED ============

int main(int argc, char *argv[]){

    // ========================= READING INPUT VARIABLES ==========================

    if (argc <= 5){
        printf("Missing inputs!!! Need ATTACK ANGLE INDEX, WIND X, WIND Y, WIND Z and 1 IF TRAJECTORY FILE IS NEEDED\n");
        return 0;
    }

    int alpha_index = atoi( *(argv + 1) );

    double W[3];
    W[0] = atof( *(argv + 2) );
    W[1] = atof( *(argv + 3) );
    W[2] = atof( *(argv + 4) );

    bool file_needed = atof( *(argv + 4) );

    if (alpha_index >= n_alphas){
        printf("Alpha index too big!!!\n");
        return 0;
    }

     // ========================= CREATING TRAJECTORY OUTPUT FILE ==========================

    char text[30];
    time_t now = time(NULL);
    struct tm tim;
    tim = *(localtime(&now));
    strftime(text, sizeof(text)-1, "%b-%d-%Y_%H-%M-%S", &tim);
    text[30] = 0;

    // concat the date to file name
    char *filename_trajectory;
    if((filename_trajectory = malloc(strlen("filename.txt")+strlen(text)+1)) != NULL){
        filename_trajectory[0] = '\0';   // ensures the memory is an empty string
        strcat(filename_trajectory,"trajectory-");
        strcat(filename_trajectory,text);
        strcat(filename_trajectory,".txt");
    }

    FILE *trajectory, *wind;
    trajectory = fopen("out.txt", "w+"); // fopen(filename_trajectory, "w+");

    fprintf(trajectory, "t       x_kite       y_kite       z_kite      x_blocco      y_blocco       z_blocco      theta      phi      r_diff\n");

    // ============================ VARIABLES DEFINITION ============================

    // kite motion vectors from fixed origin (x, z)
    double *rk = (double*) malloc(dim * sizeof(double)); 
    double *vk = (double*) malloc(dim * sizeof(double)); 
    double *ak = (double*) malloc(dim * sizeof(double)); 

    // block motion vectors from fixed origin  (x, z)
    double *r_block = (double*) malloc(dim * sizeof(double)); 
    double *v_block = (double*) malloc(dim * sizeof(double));
    double *a_block = (double*) malloc(dim * sizeof(double));  

    // block motion vectors from fixed origin  (x, z)
    double *r_diff = (double*) malloc(dim * sizeof(double)); 
    double *v_diff = (double*) malloc(dim * sizeof(double));
    double *a_diff = (double*) malloc(dim * sizeof(double)); 

    double theta = _theta0;
    double phi = phi0;
    double dtheta = _dtheta0;
    double dphi = dphi0;
    double r_diff_modulo;
    double v_diff_modulo;

    double lift=0, drag=0;
    double T = 0;
    double F_attr;
    double F_vinc;
    double theta_star;
    int stability = 0;
    int decollato = 0;

    int sector = 0;

    variables_initialization(rk, vk, ak, theta, phi, dtheta, dphi, r_block, v_block, a_block, r_diff, v_diff, a_diff);

    printf("init: %f       %f      %f      %f      %f      %f\n\n", rk[0], rk[1], rk[2], r_block[0], r_block[1], r_block[2]);

    int t = 0;   

    // ============================ KITE FLYING LOOP, STOPS WHEN FALL OCCURS ============================
    
    for (int i=0; i<STEPS; i++){

        integration_trajectory(rk, vk, ak, r_block, v_block, a_block, r_diff, v_diff, a_diff, \
                            &theta, &phi, alpha_index, mu, W, &lift, &drag, &T, &F_attr, i, &sector);

        r_diff_modulo = sqrt(r_diff[0]*r_diff[0] + r_diff[1]*r_diff[1] + r_diff[2]*r_diff[2]);
        F_vinc = m_block*g - T*sin(theta);

        // moving the kite to put it again at distance R with the block

        rk[0] = r_block[0] + (rk[0] - r_block[0])/fabs(r_diff_modulo)*R;
        rk[1] = r_block[1] + (rk[1] - r_block[1])/fabs(r_diff_modulo)*R;
        rk[2] = r_block[2] + (rk[2] - r_block[2])/fabs(r_diff_modulo)*R;

        if (i%PRINTSTEP == 0 || rk[2] <= 0.){
            fprintf(trajectory, "%d       %f       %f      %f      %f      %f      %f      %f      %f      %f\n", \
                    t, rk[0], rk[1], rk[2], r_block[0], r_block[1], r_block[2], theta, phi, r_diff_modulo);
            if (rk[2] <=0. ){ // maybe put rk[2] = 0
                printf("Kite Fall, steps %d, z<0, break\n", i);
                break;
            }
        }

        t += 1;

        if (F_vinc < 0) {
            decollato = 1;
        }

    }

    v_diff_modulo = sqrt(v_diff[0]*v_diff[0] + v_diff[1]*v_diff[1] + v_diff[2]*v_diff[2]);

    dtheta = -1/(sqrt(1-(r_diff[2]/r_diff_modulo)*(r_diff[2]/r_diff_modulo)))*(v_diff[2]*r_diff_modulo-r_diff[2]*v_diff_modulo)/(r_diff_modulo*r_diff_modulo);

    dphi = 1/(1+(r_diff[1]/r_diff[0])*(r_diff[1]/r_diff[0]))*(v_diff[1]*r_diff[0] - r_diff[0]*v_diff[1])/(r_diff[0]*r_diff[0]);
    
    if (rk[2] <= 0){
        rk[2] = 0;
        v_block[0] = 0;
        v_block[1] = 0;
        vk[0] = 0;
        vk[1] = 0;
        vk[2] = 0;
        theta = PI/2;
        dtheta = 0;
    }

    printf("iter, tot time, m_block, alpha, mu, theta0, theta_fin, v_block_fin_x, v_block_fin_y, Wind_x, Wind_y, Wind_z");
    printf(" vrelkite_x, vrelkite_y, vrelkite_z, F_vinc, Tension, Lift, Drag\n");
    
    printf("%d, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f\n", \
    t, t*h, m_block, alphas[alpha_index], mu, _theta0, theta, v_block[0], v_block[1], W[0], W[1], W[2], \
    vk[0] - W[0], vk[1] - W[1], vk[2] - W[2], F_vinc, T, lift, drag);

    free(rk);
    free(vk);
    free(ak);

    free(r_block);
    free(v_block);
    free(a_block);

    free(r_diff);
    free(v_diff);
    free(a_diff);

    fclose(trajectory);
    if (file_needed == 0){
        remove(filename_trajectory);
    }
    //remove("a.out");

    return 0;

    }

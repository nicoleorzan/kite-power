#include "Dynamics/dynamics_3d_cartesian.h"
//#include "Dynamics/winds.h"
#include <time.h>
#include <string.h>
#include <stdbool.h>

#define STEPS 1
#define n_alphas 16
#define theta0 PI/4
#define vtheta0 0.
#define dim 3

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
    double *rk1 = (double*) malloc(dim * sizeof(double)); 
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

    double theta = PI/4.;
    double phi = 0;
    double r_diff_modulo = sqrt(r_diff[0]*r_diff[0] + r_diff[1]*r_diff[1] + r_diff[2]*r_diff[2]);

    double lift=0, drag=0;
    double mu = -0.0872665;
    double F_vinc;
    double theta_star;
    double T = 0;
    int stability = 0;
    int decollato = 0;

    variables_initialization(rk, vk, ak, theta, phi, r_block, v_block, a_block);
    printf("init: %f       %f      %f      %f      %f      %f\n\n", rk[0], rk[1], rk[2], r_block[0], r_block[1], r_block[2]);

    int t = 0;   

    // ============================ KITE FLYING LOOP, STOPS WHEN FALL OCCURS ============================
    
    for (int i=0; i<STEPS; i++){

        integration_trajectory(rk, vk, ak, r_block, v_block, a_block, \
                            r_diff, v_diff, a_diff, &theta, &phi, alpha_index, mu, W, &lift, &drag, &T, i);

        if (m_block*g < T*cos(theta)){
            printf("m_block*g < T*cos(theta), exiting\n");
            break;
        }

        r_diff_modulo = sqrt(r_diff[0]*r_diff[0] + r_diff[1]*r_diff[1] + r_diff[2]*r_diff[2]);
        
        if (rk[2] <= 0.) {
            printf("Kite Fall, steps %d, z<0, break\n", i);
            fprintf(trajectory, "%d       %f       %f      %f      %f      %f      %f      %f      %f      %f\n", \
                    t, rk[0], rk[1], 0., r_block[0], r_block[1], r_block[2], theta, phi, r_diff_modulo);
            break;
        }

        if (i%1000 == 0){

            rk1[0] = r_block[0] + (rk[0] - r_block[0])/fabs(r_diff_modulo)*R;
            rk1[1] = r_block[1] + (rk[1] - r_block[1])/fabs(r_diff_modulo)*R;
            rk1[2] = r_block[2] + (rk[2] - r_block[2])/fabs(r_diff_modulo)*R;

            /*printf("(rk[0] - r_block[0])/fabs(r_diff_modulo)*R=%f\n", (rk[0] - r_block[0])/fabs(r_diff_modulo)*R);
            printf("(rk[1] - r_block[1])/fabs(r_diff_modulo)*R=%f\n", (rk[1] - r_block[1])/fabs(r_diff_modulo)*R);
            printf("(rk[2] - r_block[2])/fabs(r_diff_modulo)*R=%f\n", (rk[2] - r_block[2])/fabs(r_diff_modulo)*R);*/

            rk[0] = rk1[0];
            rk[1] = rk1[1];
            rk[2] = rk1[2];

            //printf("%f       %f      %f      %f      %f      %f\n", rk[0], rk[1], 0., r_block[0], r_block[1], r_block[2]);
            fprintf(trajectory, "%d       %f       %f      %f      %f      %f      %f      %f      %f      %f\n", \
                    t, rk[0], rk[1], rk[2], r_block[0], r_block[1], r_block[2], theta, phi, r_diff_modulo);
        }

        t += 1;

        F_vinc = m_block*g - T*cos(theta);

        if (F_vinc < 0) {
            decollato = 1;
        }
    }

    if ( rk[1] > 0.) {
        printf("iter, alpha, mu, theta0, Theta_fin, v_block_fin_x, v_block_fin_y, F_vinc, ");
        printf("Tension, Lift, Drag, Wind_x, Wind_y, Wind_z\n");
        
        printf("%d, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f\n", \
        t, alphas[alpha_index], mu, theta0, theta, v_block[0], v_block[1], \
        F_vinc, T, lift, drag, W[0], W[1], W[2]);
    } else {
        printf("iter, alpha, mu, theta0, Theta_fin, v_block_fin_x, v_block_fin_y, F_vinc, ");
        printf("Tension, Lift, Drag, Wind_x, Wind_y, Wind_z\n");

        printf("%d, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f\n", \
        t, alphas[alpha_index], mu, theta0, theta, v_block[0], v_block[1], \
        F_vinc, T, lift, drag, W[0], W[1], W[2]);
    }

    free(rk);
    free(rk1);
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

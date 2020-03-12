#include "Dynamics/dynamics_2d_spherical.h"
#include "Dynamics/winds.h"
#include <time.h>
#include <string.h>
#include <stdbool.h>

#define STEPS 1000000
#define theta0 0.
#define vtheta0 0.5

// ============== FILE INPUT: ATTACK ANGLE AND WIND COEFF AND IF TRAJECTORY IS NEEDED ============

int main(int argc, char *argv[]){

    // ========================= READING INPUT VARIABLES ==========================

    if (argc <= 4){
        printf("Missing inputs!!! Need ATTACK ANGLE INDEX, WIND X, WIND Y and 1 IF TRAJECTORY FILE IS NEEDED\n");
        return 0;
    }

    int alpha_index = atoi( *(argv + 1) );

    double W[2];
    W[0] = atof( *(argv + 2) );
    W[1] = atof( *(argv + 3) );

    //double wind_coeff = 0;

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

    fprintf(trajectory, "t         x_kite          z_kite         x_block         z_block       wind_x      wind_y      v_block\n");

    // ============================ VARIABLES DEFINITION ============================

    // kite motion vectors from fixed origin (x, z)
    double *rk = (double*) malloc(2 * sizeof(double)); 
    double *vk = (double*) malloc(2 * sizeof(double)); 
    double *ak = (double*) malloc(2 * sizeof(double)); 

    // block motion vectors from fixed origin  (x, z)
    double *r_block = (double*) malloc(2 * sizeof(double)); 
    double *v_block = (double*) malloc(2 * sizeof(double));
    double *a_block = (double*) malloc(2 * sizeof(double));  

    // theta, dtheta, ddtheta of kite from block sdr
    double *theta = (double*) malloc(3 * sizeof(double)); 
     
    double lift=0, drag=0;
    double T = 0;
    double F_vinc;
    double theta_star;
    int stability = 0;
    int decollato = 0;

    variables_initialization(rk, vk, ak, theta0, vtheta0, r_block, v_block, a_block, theta);

    int t = 0;   

    // ============================ KITE FLYING LOOP, STOPS WHEN FALL OCCURS ============================
    
    for (int i=0; i<STEPS; i++){

        //streamfunction(rk, &W[0], &W[1]);

        integration_trajectory(rk, vk, ak, r_block, v_block, a_block, theta, &T, alpha_index, \
                             W, &lift, &drag);
        
        if (rk[1] <= 0.) {
            printf("Kite Fall, steps %d, z<0, break\n", i);
            fprintf(trajectory, "%d       %f       %f      %f      %f      %f      %f     %f\n", \
                    t, rk[0], 0., r_block[0], r_block[1], W[0], W[1], v_block[0]);
            break;
        }
        
        theta_star = atan((lift - m*g)/drag);

        if (i%1000 == 0){
            fprintf(trajectory, "%d       %f       %f      %f      %f      %f      %f     %f\n", \
                    t, rk[0], rk[1], r_block[0], r_block[1], W[0], W[1], v_block[0]);
        }

        t += 1;

        F_vinc = m_block*g - T*sin(theta[0]);

        if (F_vinc < 0) {
            decollato = 1;
        }
    }

    stability = 0; // 0 = stability not reached

    if (theta[0] != 0. && abs(theta[0] == theta_star) < 10E-5 && theta[1] < 10E-5){
        stability = 1;
    }

    if ( rk[1] > 0.) {
        
        printf("Alpha, theta0, vtheta0, Wx, Wy, Theta_fin, Vel_theta_fin, V_block_fin, Vk-W, F_vinc, Tension, Lift, Drag, Stab, Decol\n");
        
        printf("%.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %d, %d\n", \
        alphas[alpha_index], theta0, vtheta0, W[0], W[1], theta[0], theta[1], v_block[0], \
        sqrt((vk[0] - W[0])*(vk[0] - W[0]) + (vk[1] - W[1])*(vk[1] - W[1])), \
        F_vinc, T, lift, drag, stability, decollato);
    } else {
        printf("Alpha, theta0, vtheta0, WindC, Wx, Wy, Theta_fin, Vel_theta_fin, V_block_fin,  Vk-W, F_vinc, Tension, Lift, Drag, Stab, Decol\n");
        
        printf("%.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %d, %d\n", \
        alphas[alpha_index], theta0, vtheta0, W[0], W[1], 0.00, 0.00, 0.00, \
        sqrt((0. - W[0])*(0. - W[0]) + (0. - W[1])*(0. - W[1])) , \
        F_vinc, T, lift, drag, stability, decollato);
    }

    free(rk);
    free(vk);
    free(ak);

    free(r_block);
    free(v_block);
    free(a_block);

    free(theta);

    fclose(trajectory);
    if (file_needed == 0){
        remove(filename_trajectory);
    }

    return 0;

    }

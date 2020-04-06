#include "Dynamics/dynamics_2d_polar.h"
#include "Dynamics/winds.h"
#include <time.h>
#include <string.h>
#include <stdbool.h>

#define dim 2

// ============== FILE INPUT: ATTACK ANGLE AND WIND X, Y ============

int main(int argc, char *argv[]){

    // ========================= READING INPUT VARIABLES ==========================

    if (argc <= 3){
        printf("Missing inputs!!! Need ATTACK ANGLE INDEX, WIND X, WIND Y\n");
        return 0;
    }

    int alpha_index = atoi( *(argv + 1) );

    double W[2];
    W[0] = atof( *(argv + 2) );
    W[1] = atof( *(argv + 3) );

    if (alpha_index >= n_alphas){
        printf("Alpha index too big!!!\n");
        return 0;
    }

    // ========================= CREATING TRAJECTORY OUTPUT FILE ==========================

    /*char text[30];
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
    }*/

    FILE *trajectory;
    trajectory = fopen("out.txt", "w+"); 

    fprintf(trajectory, "t,x_kite,z_kite,x_block,z_block,theta,vtheta,windx,windy,v_block,Tension\n");

    // ============================ VARIABLES DEFINITION ============================

    // kite motion vectors from fixed origin (x, z)
    double *rk = (double*) malloc(dim * sizeof(double)); 
    double *vk = (double*) malloc(dim * sizeof(double)); 
    double *ak = (double*) malloc(dim * sizeof(double)); 

    // block motion vectors from fixed origin  (x, z)
    double *r_block = (double*) malloc(dim * sizeof(double)); 
    double *v_block = (double*) malloc(dim * sizeof(double));
    double *a_block = (double*) malloc(dim * sizeof(double));  

    // theta, dtheta, ddtheta of kite from block sdr
    double *theta = (double*) malloc(3 * sizeof(double)); 
     
    double lift = 0;
    double drag = 0;
    double T = 0;
    int sector = 0;
    double F_vinc = 0;
    double F_attr = 0;
    double theta_star = 0;
    int stability = 0;
    int decollato = 0;
    
    variables_initialization(rk, vk, ak, theta0, vtheta0, r_block, v_block, a_block, theta);

    //streamfunction2d(rk, W);

    int t = 0;   

    // ============================ KITE FLYING LOOP, STOPS WHEN FALL OCCURS ============================
    
    for (int i=0; i<STEPS; i++){

        integration_trajectory(rk, vk, ak, r_block, v_block, a_block, theta, alpha_index, \
                             W, &lift, &drag, &T, &F_attr, i, &sector);

        //streamfunction2d(rk, W);  

        if (i%PRINTSTEP == 0  || rk[1] <= 0.){
            fprintf(trajectory, "%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", \
                    t, rk[0], rk[1], r_block[0], r_block[1], theta[0], theta[1], W[0], W[1], v_block[0], T);
            if (rk[1] <=0. ){
                //printf("Kite Fall, steps %d, z<0, break\n", i);
                break;
            }
        }

        F_vinc = m_block*g - T*sin(theta[0]);

        t += 1;

        if (F_vinc < 0) {
            decollato = 1;
        }
    }

    theta_star = atan((lift - m*g)/drag);

    if (theta[0] != 0. && fabs(theta[0] - theta_star) < THRESHOLD && theta[1] < THRESHOLD){
        stability = 1;
    }

    if (rk[1] <= 0){
        rk[1] = 0;
        v_block[0] = 0;
        vk[0] = 0;
        vk[1] = 0;
        theta[0] = 0;
        theta[1] = 0;
    }

    //printf("iter, tot time, m_block, alpha, theta0, theta_fin, v_theta_fin, v_block_fin_x, Wind_x, Wind_y, ");
    //printf(" vkitex, vkitey, vrelkite_x, vrelkite_y, F_vinc, Tension, Lift, Drag, Stability\n");
    
    printf("%d, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %d\n", \
    t, t*h, m_block, alphas[alpha_index], theta0, theta[0], theta[1], v_block[0], W[0], W[1], \
    vk[0], vk[1], vk[0] - W[0], vk[1] - W[1], F_vinc, T, lift, drag, stability);

    free(rk);
    free(vk);
    free(ak);

    free(r_block);
    free(v_block);
    free(a_block);

    free(theta);

    fclose(trajectory);

    return 0;

    }

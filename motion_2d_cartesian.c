#include "Dynamics/dynamics_2d_cartesian_tension.h"
//#include "Dynamics/winds.h"
#include <time.h>
#include <string.h>
#include <stdbool.h>

#define theta0 0.
#define vtheta0 0.5
#define dim 2

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

    fprintf(trajectory, "t       x_kite       z_kite      x_block      z_block      theta      r_diff      T\n");

    // ============================ VARIABLES DEFINITION ============================

    // kite motion vectors from fixed origin (x, z)
    double *rk = (double*) malloc(dim * sizeof(double)); 
    double *vk = (double*) malloc(dim * sizeof(double)); 
    double *ak = (double*) malloc(dim * sizeof(double)); 

    // block motion vectors from fixed origin  (x, z)
    double *r_block = (double*) malloc(dim * sizeof(double)); 
    double *v_block = (double*) malloc(dim * sizeof(double));
    double *a_block = (double*) malloc(dim * sizeof(double));  

    double *r_diff = (double*) malloc(dim * sizeof(double)); 
    double *v_diff = (double*) malloc(dim * sizeof(double)); 
    double *a_diff = (double*) malloc(dim * sizeof(double)); 

    double theta = theta0;

    double r_diff_modulo;

    double lift=0, drag=0;
    double T = 0;
    int sector = 0;
    double F_vinc;
    double theta_star;
    int stability = 0;
    int decollato = 0;

    variables_initialization(rk, vk, ak, theta0, vtheta0, r_block, v_block, a_block, r_diff, v_diff, a_diff);

    int t = 0;   

    // ============================ KITE FLYING LOOP, STOPS WHEN FALL OCCURS ============================
    
    for (int i=0; i<STEPS; i++){

        integration_trajectory(rk, vk, ak, r_block, v_block, a_block, r_diff, v_diff, a_diff, \
                            &theta, alpha_index, W, &lift, &drag, &T, i, &sector);
        printf("T=%f, sector=%d\n", T, sector);

        /*if (m_block*g < T*sin(theta)){
            printf("m_block*g < T*sin(theta), exiting\n");
            break;
        }*/

        r_diff_modulo = sqrt(r_diff[0]*r_diff[0] + r_diff[1]*r_diff[1]);
        theta_star = atan((lift - m*g)/drag); // e` giusto calcolarlo ad ogni step??
        F_vinc = m_block*g - T*sin(theta);

        // moving the kite to put it again at distance R with the block

        rk[0] = r_block[0] + (rk[0] - r_block[0])/fabs(r_diff_modulo)*R;
        rk[1] = r_block[1] + (rk[1] - r_block[1])/fabs(r_diff_modulo)*R;

        if (i%PRINTSTEP == 0 || rk[1] <= 0.){
            fprintf(trajectory, "%d       %f       %f      %f      %f      %f      %f      %f\n", \
                    t, rk[0], rk[1], r_block[0], r_block[1], theta, r_diff_modulo, T);
            if (rk[1] <=0. ){
                printf("Kite Fall, steps %d, z<0, break\n", i);
                break;
            }
        }

        t += 1;

        if (F_vinc < 0) {
            decollato = 1;
        }

    }

    //stability = 0; // 0 = stability not reached

    /*if (theta[0] != 0. && abs(theta[0] == theta_star) < 10E-5 && theta[1] < 10E-5){
        stability = 1;
    }*/

    printf("iter, alpha, theta0, Theta_fin, v_block_fin_x, F_vinc, ");
    printf("Tension, Lift, Drag, Wind_x, Wind_y\n");
    
    printf("%d, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f\n", \
    t, alphas[alpha_index], theta0, theta, v_block[0], \
    F_vinc, T, lift, drag, W[0], W[1]);

    free(rk);
    free(vk);
    free(ak);

    free(r_diff);
    free(v_diff);
    free(a_diff);

    free(r_block);
    free(v_block);
    free(a_block);

    fclose(trajectory);
    if (file_needed == 0){
        remove(filename_trajectory);
    }
    //remove("a.out");

    return 0;

}

#include "Dynamics/dynamics_2d_cartesian.h"
//#include "Dynamics/winds.h"
#include <time.h>
#include <string.h>
#include <stdbool.h>

#define dim 2

// ============== FILE INPUT: ATTACK ANGLE AND WIND X, Z ============

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

    FILE *trajectory, *debug;//, *debug_et;
    trajectory = fopen("outc.txt", "w+");
    debug = fopen("debug2d.csv", "w+");
    // debug_et = fopen("Et-analysis/debug_et.csv", "w+");

    fprintf(trajectory, "t,x_kite,z_kite,x_block,z_block,theta,vtheta,windx,windy,v_block,Tension\n");
    fprintf(debug, "i,Alpha,theta,Windx,Windz,Vkx,Vkz,vk[0]-W[0],vk[1]-W[1],Lift,Liftx,Liftz,Drag,Tension,F_attrito,sector\n");
    //fprintf(debug_et, "i,Alpha,Windx,Windz,vkx,vkz,t2_value\n");

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
    double vtheta = vtheta0;

    double r_diff_modulo = 0;

    double lift = 0;
    double drag = 0;
    double T = 0;
    int sector = 0;
    double F_vinc = 0;
    double F_attr = 0;
    double theta_star = 0;
    int stability = 0;
    int decollato = 0;
    double dtheta = 0;

    double l0, l1;
    double d0, d1;

    variables_initialization(rk, vk, ak, theta, vtheta, r_block, v_block, a_block, r_diff, v_diff, a_diff);
    
    //streamfunction2d(rk, W);

    int t = 0;   
    int et_val = 0;

    // ============================ KITE FLYING LOOP, STOPS WHEN FALL OCCURS ============================
    
    for (int i=0; i<STEPS; i++){

        integration_trajectory(rk, vk, ak, r_block, v_block, a_block, r_diff, v_diff, a_diff, \
                            &theta, &vtheta, alpha_index, W, &lift, &drag, &T, &F_attr, i, &sector, &l0, &l1, &d0, &d1, &et_val);
        r_diff_modulo = sqrt(r_diff[0]*r_diff[0] + r_diff[1]*r_diff[1]);
        theta_star = atan((lift - m*g)/drag);
        F_vinc = m_block*g - T*sin(theta);

        // moving the kite to put it again at distance R with the block

        rk[0] = r_block[0] + (rk[0] - r_block[0])/fabs(r_diff_modulo)*R;
        rk[1] = r_block[1] + (rk[1] - r_block[1])/fabs(r_diff_modulo)*R;

        //fprintf(debug_et, "%d,%.2f,%.2f,%.2f,%.2f,%.2f,%d\n", i, alphas[alpha_index], W[0], W[1], vk[0], vk[1], et_val);

        //streamfunction2d(rk, W);

        if (i%1 == 0){
            fprintf(debug, "%d,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%d\n",
                    i, alphas[alpha_index], theta, W[0], W[1], vk[0], vk[1], vk[0]-W[0], vk[1]-W[1], lift, l0, \
                    l1, drag, T, fabs(F_attr), sector);
        }

        if (i%PRINTSTEP == 0 || rk[1] <= 0.){
            //fprintf(trajectory, "%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", \
            //         t, rk[0], rk[1], r_block[0], r_block[1], theta, vtheta, W[0], W[1], v_block[0], T);
            if (rk[1] <=0. ){
                //printf("Kite Fall, steps %d, z<0, break\n", i);
                break;
            }
        }

        t += 1;

        if (F_vinc < 0) {
            decollato = 1;
        }
    }

    theta_star = atan((lift - m*g)/drag);

    if (theta != 0. && abs(theta - theta_star) < THRESHOLD && dtheta < THRESHOLD){
        stability = 1;
    }

    if (rk[1] <= 0){
        rk[1] = 0;
        v_block[0] = 0;
        vk[0] = 0;
        vk[1] = 0;
        theta = 0;
        dtheta = 0;
    }

    //printf("\niter, tot time, m_block, alpha, theta0, theta_fin, v_theta_fin, v_block_fin_x, Wind_x, Wind_y, ");
    //printf(" vkitex, vkitez, vrelkite_x, vrelkite_y, F_vinc, Tension, Lift, Drag, Stability\n\n");
    
    printf("%d, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %d\n", \
    t, t*h, m_block, alphas[alpha_index], theta0, theta, dtheta, v_block[0], W[0], W[1], \
    vk[0], vk[1], vk[0] - W[0], vk[1] - W[1], F_vinc, T, lift,  drag, stability);

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
    fclose(debug);
    //fclose(debug_et);

    return 0;

}

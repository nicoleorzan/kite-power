#include "../../Dynamics/dynamics_2d_cartesian.h"
#include "../../Dynamics/winds.h"
#include "sarsa_alpha.h"

double reward_dt = h*decision_time;
#define dim 2

int main(int argc, char *argv[]){

    FILE *out;
    out = fopen("out_eval.txt", "w");

    fprintf(out, "step,x_kite,z_kite,x_block,z_block,theta,vtheta,windx,windy,v_block,vkx,vkz,alpha,action,reward,tot_return,lx,lz,dx,dz,lxdx,lzdz,vkx,vkz,vkmod\n");
    // ======== DYNAMICS VARIABLES =======

    // vettori moto kite dall'origine fissa (x, z)
    double *rk = (double*) malloc(dim * sizeof(double)); 
    double *vk = (double*) malloc(dim * sizeof(double)); 
    double *ak = (double*) malloc(dim * sizeof(double)); 

    // vettori moto block dall'origine fissa (x, z)
    double *r_block = (double*) malloc(dim * sizeof(double));
    double *v_block = (double*) malloc(dim * sizeof(double));
    double *a_block = (double*) malloc(dim * sizeof(double));

    double *r_diff = (double*) malloc(dim * sizeof(double));
    double *v_diff = (double*) malloc(dim * sizeof(double));
    double *a_diff = (double*) malloc(dim * sizeof(double));

    double r_diff_modulo;
    
    double theta;
    double vtheta;

    double W[dim] = {10, 0};
    int et_val = 0;

    double T = 0;
    double F_attr = 0;
    double lift=0, drag=0;
    int sector = 0;
    double l0, l1;
    double d0, d1;

    // ======== LEARNING VARIABLES =======

    double * Q = (double*) malloc(n_alphas * n_actions * sizeof(double));

    double reward = 0;
    double tot_reward = 0;

    // SARSA STATE
    int s_alpha;

    // SARSA ACTION
    int a_alpha = 0;           

    int it = 0;

    rk[0] = 0;
    rk[1] = R;

    streamfunction2d_hard(rk, W);
    printf("W0 = %f, W1 = %f\n", W[0], W[1]);
    double max_return = sqrt(W[0]*W[0] + W[1]*W[1])*max_steps*h;
    printf("reward max %f\n\n", max_return);


    load_matrix(Q, "matrix_streamfunction.txt");
    //print_mat(Q);
    
    // ======================= EPISODE INITIALIZATION ==========================

    //printf("Cartesian code\n");

    theta = PI/4; //theta0;
    vtheta = 0.; //vtheta0;

    variables_initialization(rk, vk, ak, theta, vtheta, r_block, v_block, a_block, r_diff, v_diff, a_diff);

    streamfunction2d_hard(rk, W);

    s_alpha = 12; //s_alpha0;
    a_alpha = 1; //select_action_greedy(Q, s_alpha);

    it = 0;
    reward = 0;
    tot_reward = 0.; 

    printf("theta0 = %f, vel_theta0 = %f\n", theta, vtheta);
    printf("initial alpha state=%d, initial action=%d\n", s_alpha, a_alpha);
    printf("W[0]=%f, W[1]=%f\n", W[0], W[1]);

    integration_trajectory(rk, vk, ak, r_block, v_block, a_block, r_diff, v_diff, a_diff, \
                        &theta, &vtheta, s_alpha, W, &lift, &drag, &T, &F_attr, it, &sector, &l0, &l1, &d0, &d1, &et_val);

    streamfunction2d_hard(rk, W);

    while (rk[1] > 0){

        if (it%decision_time == 0){
            fprintf(out, "%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%fx\n", it, rk[0], rk[1], 
                r_block[0], r_block[1], theta, vtheta, W[0], W[1], v_block[0], vk[0], vk[1], alphas[s_alpha], a_alpha, reward,
                tot_reward, l0, l1, d0, d1, l0+d0, l1+d1, vk[0], vk[1], sqrt(vk[0]*vk[0] + vk[1]*vk[1]));
            //printf("vk[0]=%f, vk[1]=%f\n", vk[0], vk[1]);
            //printf("dist = %f\n", sqrt((rk[0]-r_block[0])*(rk[0]-r_block[0]) +  (rk[1]-r_block[1])*(rk[1]-r_block[1])));
        }   

        if (it > max_steps){
            printf("MAX STEPS, %d, exiting\n", max_steps);
            printf("final alpha=%d\n", s_alpha);
            printf("return=%f, space percurred=%f\n\n", tot_reward, r_block[0]);

            break;
        }

        integration_trajectory(rk, vk, ak, r_block, v_block, a_block, r_diff, v_diff, a_diff, \
                        &theta, &vtheta, s_alpha, W, &lift, &drag, &T, &F_attr, it, &sector, &l0, &l1, &d0, &d1, &et_val);

        r_diff_modulo = sqrt(r_diff[0]*r_diff[0] + r_diff[1]*r_diff[1]);

        rk[0] = r_block[0] + (rk[0] - r_block[0])/fabs(r_diff_modulo)*R;
        rk[1] = r_block[1] + (rk[1] - r_block[1])/fabs(r_diff_modulo)*R;

        streamfunction2d_hard(rk, W);

        reward = fabs(v_block[0])*reward_dt;

        if (it%decision_time == 0){
            tot_reward += reward;
        }
        
        if (rk[1] <= 0.) {

            printf("Kite fallen: z<0, steps=%d, break\n", it);
            printf("final alpha=%d\n", s_alpha);
            printf("return-penalty=%f, space percurred=%f\n\n", tot_reward + PENALTY, r_block[0]);

            fprintf(out, "%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%fx\n", it, rk[0], rk[1], 
                r_block[0], r_block[1], theta, vtheta, W[0], W[1], v_block[0], vk[0], vk[1], alphas[s_alpha], a_alpha, reward,
                tot_reward, l0, l1, d0, d1, l0+d0, l1+d1, vk[0], vk[1], sqrt(vk[0]*vk[0] + vk[1]*vk[1]));
            break;
        }

        if (it%decision_time == 0){

            // FIND STATE S1 USING ACTION A
            s_alpha = update_state(s_alpha, a_alpha);

            // SEARCH FOR NEXT ACTION A1
            a_alpha = select_action_greedy(Q, s_alpha);
            
        }

        it += 1;
    }

    free(rk);
    free(vk);
    free(ak);

    free(r_diff);
    free(v_diff);
    free(a_diff);

    free(r_block);
    free(v_block);
    free(a_block);

    free(Q);

    fclose(out);

    return 0;

    }

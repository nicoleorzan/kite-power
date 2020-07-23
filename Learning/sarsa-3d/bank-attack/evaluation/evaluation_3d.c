#include "../../../../Dynamics/dynamics_3d_cartesian_rail.h"
#include "evaluation_3d.h"

double reward_dt = h*decision_time;

#define _theta0 PI/4.
#define _phi0 0.
#define _dtheta0 0.
#define _dphi0 0.
#define dim 3

int main(int argc, char *argv[]){

    FILE *out;
    out = fopen("out_power.txt", "w");

    fprintf(out, "step,x_kite,y_kite,z_kite,x_block,y_block,z_block,theta,phi,windx,windy,windz,vrelx,vrely,vrelz,v_blockx,v_blocky,alpha,alpha_action,mu,mu_action,reward,l0,l1,l2,d0,d1,d2,f_attr\n");
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
    double phi;
    double vphi;

    double T = 0;
    double F_attr = 0;
    double lift = 0;
    double drag = 0;
    int sector = 0;
    double l0, l1, l2;
    double d0, d1, d2;

    // ======== LEARNING VARIABLES =======

    double * Q = (double*) malloc(n_alphas * n_bank * n_actions * n_actions * sizeof(double));

    double reward = 0;
    double tot_reward = 0;

    // SARSA STATE
    int s_alpha, s_alpha1;
    int s_bank, s_bank1;

    // SARSA ACTION
    int a_alpha = 0, a_alpha1 = 0;  
    int a_bank = 0, a_bank1 = 0;          

    int it = 0;

    double W[dim] = {10, 0, 0};

    printf("load matrix\n");
    load_matrix_3d(Q, "final_matrix.txt");

    print_mat(Q);
    // ======================= EPISODE INITIALIZATION ==========================

    printf("Cartesian code\n");

    theta = _theta0;
    vtheta = _dtheta0;
    phi = _phi0;
    vphi = _dphi0;

    variables_initialization(rk, vk, ak, theta, phi, vtheta, vphi, r_block, v_block, a_block, r_diff, v_diff, a_diff);

    s_alpha = s_alpha0;
    s_bank = s_bank0;
    select_action_greedy(Q, s_alpha, &a_alpha, s_bank, &a_bank);

    it = 0;
    reward = 0;
    tot_reward = 0.; 

    printf("theta0 = %f, vel_theta0 = %f\n", theta, vtheta);
    printf("phi0 = %f, vel_phi0 = %f\n", phi, vphi);
    printf("decision time=%d\n", decision_time);
    printf("initial state: s_alpha=%d, s_bank=%d\n", s_alpha, s_bank);
    printf("with values: alpha=%f, bank=%f\n", alphas[s_alpha], bank[s_bank]);
    printf("starting actions: a_alpha=%d, a_bank=%d\n", a_alpha, a_bank);
    printf("W[0]=%f, W[1]=%f, W[2]=%f\n", W[0], W[1], W[2]);

    integration_trajectory(rk, vk, ak, r_block, v_block, a_block, r_diff, v_diff, a_diff, \
                        &theta, &phi, s_alpha, bank[s_bank], W, &lift, &drag, &T, &F_attr, it, &sector, \
                        &l0, &l1, &l2, &d0, &d1, &d2);

    while (rk[2] > 0){

        if ( it%decision_time == 0 ){

            fprintf(out, "%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%d,%f,%d,%f,%f,%f,%f,%f,%f,%f,%f\n", it, rk[0], rk[1], rk[2],
                r_block[0], r_block[1], r_block[2], theta, phi, W[0], W[1], W[2], W[0]-vk[0], W[1]-vk[1], W[2]-vk[2], v_block[0], v_block[1], alphas[s_alpha], a_alpha, bank[s_bank], a_bank, \
                reward, l0, l1, l2, d0, d1, d2, F_attr);
        }   

        if (it > max_steps){
            printf("MAX STEPS, %d, exiting\n", max_steps);
            printf("final alpha=%d, bank=%d\n", s_alpha, s_bank);
            printf("return=%f, space percurred=%f\n\n", tot_reward, sqrt(r_block[0]*r_block[0] + r_block[1]*r_block[1]));

            break;
        }

        integration_trajectory(rk, vk, ak, r_block, v_block, a_block, r_diff, v_diff, a_diff, \
                        &theta, &phi, s_alpha, bank[s_bank], W, &lift, &drag, &T, &F_attr, it, &sector, \
                        &l0, &l1, &l2, &d0, &d1, &d2);

        r_diff_modulo = sqrt(r_diff[0]*r_diff[0] + r_diff[1]*r_diff[1] + r_diff[2]*r_diff[2]);

        rk[0] = r_block[0] + (rk[0] - r_block[0])/fabs(r_diff_modulo)*R;
        rk[1] = r_block[1] + (rk[1] - r_block[1])/fabs(r_diff_modulo)*R;
        rk[2] = r_block[2] + (rk[2] - r_block[2])/fabs(r_diff_modulo)*R;

        reward = sqrt(v_block[0]*v_block[0] + v_block[1]*v_block[1])*reward_dt;

        if (it%decision_time == 0){
            tot_reward += reward;
        }
        
        // SOME CHECKS
        if (check(s_alpha, s_bank, a_alpha, a_bank) == 1) {
            printf("CHECK: s_alpha %d\n", s_alpha);
            printf("Matrix index error!!\n");
            break;
        }
        if (a_alpha > 2 || a_alpha < 0){
            printf("a_alpha:%d\n", a_alpha);
            break;
        }
        if (a_bank > 2 || a_bank < 0){
            printf("a_alpha:%d\n", a_alpha);
            printf("ERROR IN UPDATE STATE");
            break;
        }

        if (rk[2] <= 0.) {

            fprintf(out, "%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%d,%f,%d,%f,%f,%f,%f,%f,%f,%f,%f\n", it, rk[0], rk[1], rk[2],
                r_block[0], r_block[1], r_block[2], theta, phi, W[0], W[1], W[2], W[0]-vk[0], W[1]-vk[1], W[2]-vk[2], v_block[0], v_block[1], alphas[s_alpha], a_alpha, bank[s_bank], a_bank, \
                reward, l0, l1, l2, d0, d1, d2, F_attr);

            printf("Kite fallen: z<0, steps=%d, break\n", it);
            printf("final alpha=%d, bank=%d\n", s_alpha, s_bank);
            printf("return=%f, space percurred=%f\n\n", tot_reward, sqrt(r_block[0]*r_block[0] + r_block[1]*r_block[1]));

            break;
        }

        if (it%decision_time == 0){

            // FIND STATE S1 USING ACTION A

            s_alpha1 = update_state(s_alpha, a_alpha);
            s_bank1 = update_state(s_bank, a_bank);

            // SEARCH FOR NEXT ACTION A1

            select_action_greedy(Q, s_alpha1, &a_alpha1, s_bank1, &a_bank1);

            // MOVE ON: S = S1, A = A1

            s_alpha = s_alpha1;
            s_bank = s_bank1;

            a_alpha = a_alpha1;
            a_bank = a_bank1;
            
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

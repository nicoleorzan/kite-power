#include "../../Dynamics/dynamics_3d_cartesian_rail.h"
#include "3devaluation_head.h"

double reward_dt = h*decision_time;

#define _theta0 PI/4.
#define _phi0 0.
#define _dtheta0 0.
#define _dphi0 0.
#define dim 3

int main(int argc, char *argv[]){

    FILE *out;
    out = fopen("out3d_evaluation_power.txt", "w");

    fprintf(out, "step,x_kite,y_kite,z_kite,x_block,y_block,z_block,vkx,vky,vkz,windx,windy,windz,vrelx,vrely,vrelz,v_blockx,v_blocky,alpha,mu,alpha_action,mu_action,vrel_angle_idx,return_up_now,reward,theta,phi,l0,l1,l2,d0,d1,d2,f_attr\n");
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

    double vrel_x, vrel_y, vrel_z;
    double vrel_angle;

    // ======== LEARNING VARIABLES =======

    double * Q = (double*) malloc(n_alphas * n_bank * n_actions * n_actions * n_angles * sizeof(double));

    double reward = 0;
    double tot_reward = 0;

    // SARSA STATE
    int s_alpha, s_alpha1;
    int s_bank, s_bank1;
    int s_vrel_angle, s_vrel_angle1;  

    // SARSA ACTION
    int a_alpha = 0, a_alpha1 = 0;  
    int a_bank = 0, a_bank1 = 0;        

    int it = 0;
    int episode = 0;

    double W[dim] = {10, 0, 0};

    printf("W[0]=%f, W[1]=%f, W[2]=%f\n", W[0], W[1], W[2]);
    double max_return = sqrt(W[0]*W[0] + W[1]*W[1] + W[2]*W[2])*max_steps*h;
    printf("reward max %f\n\n", max_return);

    printf("load matrix\n");
    load_matrix_3d(Q, "last_matrix.txt");

    print_mat(Q);

    // ======================= EPISODE INITIALIZATION ==========================

    printf("Cartesian code\n");

    theta = _theta0;
    vtheta = _dtheta0;
    phi = _phi0;
    vphi = _dphi0;

    variables_initialization(rk, vk, ak, theta, phi, vtheta, vphi, r_block, v_block, a_block, r_diff, v_diff, a_diff);

    vrel_x = vk[0] - W[0];
    vrel_y = vk[1] - W[1];
    vrel_z = vk[2] - W[2];

    vrel_angle = atan2(vrel_z, vrel_x);

    s_vrel_angle = find_state_angle( vrel_angle );

    s_alpha = s_alpha0;
    s_bank = s_bank0;
    select_action_greedy(Q, s_alpha, &a_alpha, s_bank, &a_bank, s_vrel_angle);

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
                fprintf(out, "%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%d,%d,%d,%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", it, rk[0], rk[1], rk[2],
                r_block[0], r_block[1], r_block[2], vk[0], vk[1], vk[2], W[0], W[1], W[2], vrel_x, vrel_y, vrel_z, v_block[0], v_block[1], \
                s_alpha, s_bank, a_alpha, a_bank, s_vrel_angle, tot_reward, reward, theta, phi, l0, l1, l2, d0, d1, d2, F_attr);
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

        //streamfunction3d_hard(rk, W);

        reward = sqrt(v_block[0]*v_block[0] + v_block[1]*v_block[1])*reward_dt;

        if (it%decision_time == 0){
            tot_reward += reward;
        }
        
        // SOME CHECKS
        if (check(s_alpha, s_bank, a_alpha, a_bank, s_vrel_angle) == 1) {
            printf("CHECK: s_alpha %d\n", s_alpha);
            printf("Matrix index error!!\n");
            break;
        }
        if (a_alpha > 2 || a_alpha < 0){
            printf("a_alpha:%d\n", a_alpha);
            printf("ERROR IN UPDATE STATE");
            break;
        }
        if (a_bank > 2 || a_bank < 0){
            printf("a_alpha:%d\n", a_alpha);
            printf("ERROR IN UPDATE STATE");
            break;
        }

        if (rk[2] <= 0.) {

            printf("Kite fallen: z<0, steps=%d, break\n", it);
            printf("final alpha=%d, bank=%d\n", s_alpha, s_bank);
            printf("return=%f, space percurred=%f\n\n", tot_reward, sqrt(r_block[0]*r_block[0] + r_block[1]*r_block[1]));

            if ( it%decision_time == 0 ){
                fprintf(out, "%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%d,%d,%d,%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", it, rk[0], rk[1], rk[2],
                r_block[0], r_block[1], r_block[2], vk[0], vk[1], vk[2], W[0], W[1], W[2], vrel_x, vrel_y, vrel_z, v_block[0], v_block[1], \
                s_alpha, s_bank, a_alpha, a_bank, s_vrel_angle, tot_reward, reward, theta, phi, l0, l1, l2, d0, d1, d2, F_attr);
            }

            break;
        }

        if (it%decision_time == 0){

            // FIND STATE S1 USING ACTION A

            s_alpha1 = update_state(s_alpha, a_alpha);
            s_bank1 = update_state(s_bank, a_bank);

            vrel_x = vk[0] - W[0];
            vrel_y = vk[1] - W[1];
            vrel_z = vk[2] - W[2];

            vrel_angle = atan2(vrel_z, vrel_x);
            s_vrel_angle1 = find_state_angle(vrel_angle);

            // SEARCH FOR NEXT ACTION A1

            select_action_greedy(Q, s_alpha1, &a_alpha1, s_bank1, &a_bank1, s_vrel_angle1);

            // MOVE ON: S = S1, A = A1
            s_alpha = s_alpha1;
            s_bank = s_bank1;
            s_vrel_angle = s_vrel_angle1;

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

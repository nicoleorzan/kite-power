#include "../../../../Dynamics/dynamics_3d_cartesian_rail.h"
#include "sarsa_alpha_with_added_mu.h"

double reward_dt = h*decision_time;

#define _theta0 PI/4.
#define _phi0 0.
#define _dtheta0 0.
#define _dphi0 0.
#define dim 3

int main(int argc, char *argv[]){

    if (num_saved_matrices > learning_episodes){
        printf("===> Error: check save matrix step\n");
        return 1;
    }

    int save_matrix_step = (int)learning_episodes/num_saved_matrices;

    FILE *out ,*rew, *Q_mat, *Q_mat_count, *policy;
    out = fopen("out3dmu_newcoeff_6.txt", "w");
    rew = fopen("crewards3dmu_newcoeff_6.txt", "w");
    Q_mat = fopen("cQ_matrix3dmu_newcoeff_6.txt", "w");
    Q_mat_count = fopen("cQ_counter3dmu_newcoeff_6.txt", "w");
    policy = fopen("cpolicy3dmu_newcoeff_6.txt", "w");
    fprintf(rew, "episode,epsilon,Alpha,steps,decision_time,return\n");
    fprintf(out, "epsiode,t,x_kite,y_kite,z_kite,x_block,y_block,z_block,theta,vtheta,phi,vphi,windx,windy,windz,v_blockx,v_blocky\n");
    fprintf(Q_mat, "episode,alpha_idx,mu_idx,action_alpha,action_mu,Q_value\n");
    fprintf(Q_mat_count, "episode,alpha_idx,mu_idx,action_alpha,action_mu,Q_count_value\n");
    fprintf(policy, "episode,step,alpha,alpha_action,mu,mu_action,reward,Q_value\n");

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
    int * Q_count = (int*) malloc(n_alphas * n_bank * n_actions * n_actions * sizeof(int));

    double reward = 0;
    double tot_reward = 0;

    // SARSA STATE
    int s_alpha, s_alpha1;
    int s_bank, s_bank1;

    // SARSA ACTION
    int a_alpha = 0, a_alpha1 = 0;  
    int a_bank = 0, a_bank1 = 0;          

    int it = 0;
    int episode = 0;

    double update;

    double W[dim] = {10, 0, 0};

    printf("W[0]=%f, W[1]=%f, W[2]=%f\n", W[0], W[1], W[2]);
    double max_return = sqrt(W[0]*W[0] + W[1]*W[1] + W[2]*W[2])*max_steps*h;
    printf("reward max %f\n\n", max_return);

    initialize_Q(Q, 1000);
    printf("load matrix\n");
    //load_matrix(Q, "final_matrix.txt");

    print_mat(Q);

    initialize_Q_count(Q_count, 0);

    printf("Total number of episodes: %d\n", learning_episodes);

    // ======================= STARTING EPISODES ==========================

    while (episode < learning_episodes){

        if (episode == (int)(learning_episodes/2)){
            Alpha = Alpha*0.1;
            printf("Decreasing learning rate: %f\n", Alpha);
        }
        if (episode == (int)(learning_episodes*3/4)){
            Alpha = Alpha*0.1;
            printf("Decreasing learning rate: %f\n", Alpha);
        }
        if (episode == (int)learning_episodes - 5000){
            epsilon = 0.001;
        }
        if (episode == (int)learning_episodes - 150){
            epsilon = 0.0;
            printf("Last episode: greedy\n");
        }
    
        // ======================= EPISODE INITIALIZATION ==========================

        printf("Cartesian code\n");
        printf("Episode: %d\nepsilon: %f\n", episode, epsilon);
        printf("Learning rate: %f\n", Alpha);

        theta = _theta0;
        vtheta = _dtheta0;
        phi = _phi0;
        vphi = _dphi0;

        variables_initialization(rk, vk, ak, theta, phi, vtheta, vphi, r_block, v_block, a_block, r_diff, v_diff, a_diff);

        s_alpha = s_alpha0;
        s_bank = s_bank0;
        a_alpha = 1;
        a_bank = 1; //select_action(epsilon, Q, s_alpha, &a_alpha, s_bank, &a_bank);

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

            if ( (episode%save_matrix_step == 0 || episode == learning_episodes-1) && it%decision_time == 0){
                fprintf(policy, "%d,%d,%f,%d,%f,%d,%f,%f\n", \
                episode, it, alphas[s_alpha], a_alpha, bank[s_bank], a_bank, reward, \
                Q[s_bank*n_actions*n_actions*n_alphas + a_bank*n_actions*n_alphas + s_alpha*n_actions + a_alpha]);

                fprintf(out, "%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", episode, it, rk[0], rk[1], rk[2],
                    r_block[0], r_block[1], r_block[2], theta, vtheta, phi, vphi, W[0], W[1], W[2], v_block[0], v_block[1]);
            }   

            if (it > max_steps){
                printf("MAX STEPS, %d, exiting\n", max_steps);
                printf("final alpha=%d, bank=%d\n", s_alpha, s_bank);
                printf("return=%f, space percurred=%f\n\n", tot_reward, sqrt(r_block[0]*r_block[0] + r_block[1]*r_block[1]));

                fprintf(rew, "%d,%f,%f,%d,%d,%f\n", episode, epsilon, Alpha, it, decision_time, tot_reward);

                if ( episode%save_matrix_step == 0 || episode == learning_episodes-1 ){
                    fill_Q_mat(Q_mat, Q, episode);
                }
                if (episode == learning_episodes-1){
                    fill_Q_count(Q_mat_count, Q_count, episode);
                }

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

            /*if (rk[2] < 10){
                reward -= 10*reward_dt;
            }*/

            if (it%decision_time == 0){
                tot_reward += reward;
            }
            
            // SOME CHECKS
            if (check(s_alpha, s_bank, a_alpha, a_bank) == 1) {
                printf("CHECK: s_alpha %d\n", s_alpha);
                printf("Matrix index error!!\n");
                episode = learning_episodes;
                break;
            }
            if (a_alpha > 2 || a_alpha < 0){
                printf("a_alpha:%d\n", a_alpha);
                printf("ERROR IN UPDATE STATE");
                episode = learning_episodes;
                break;
            }
            if (a_bank > 2 || a_bank < 0){
                printf("a_alpha:%d\n", a_alpha);
                printf("ERROR IN UPDATE STATE");
                episode = learning_episodes;
                break;
            }

            if (rk[2] <= 0.) {

                fprintf(rew, "%d,%f,%f,%d,%d,%f\n", episode, epsilon, Alpha, it, decision_time, tot_reward);

                printf("Kite fallen: z<0, steps=%d, break\n", it);
                printf("final alpha=%d, bank=%d\n", s_alpha, s_bank);
                printf("return=%f, space percurred=%f\n\n", tot_reward, sqrt(r_block[0]*r_block[0] + r_block[1]*r_block[1]));

                //update = expected_update(Q);

                Q[s_bank*n_actions*n_alphas*n_actions + a_bank*n_actions*n_alphas + s_alpha*n_actions + a_alpha] +=
                    Alpha*(reward + PENALTY - \
                    Q[s_bank*n_actions*n_actions*n_alphas + a_bank*n_actions*n_alphas +  s_alpha*n_actions + a_alpha]);

                Q_count[s_bank*n_actions*n_alphas*n_actions + a_bank*n_actions*n_alphas + s_alpha*n_actions + a_alpha] += 1;

                if ( episode%save_matrix_step == 0 || episode == learning_episodes-1 ){
                    fill_Q_mat(Q_mat, Q, episode);
                }
                if (episode == learning_episodes-1){
                    fill_Q_count(Q_mat_count, Q_count, episode);
                }

                break;
            }

            if (it%decision_time == 0){

                // FIND STATE S1 USING ACTION A

                s_alpha1 = update_state(s_alpha, a_alpha);
                s_bank1 = update_state(s_bank, a_bank);

                // SEARCH FOR NEXT ACTION A1

                select_action(epsilon, Q, s_alpha1, &a_alpha1, s_bank1, &a_bank1);

                /*printf("%f, %f, %f\n", Q[s_bank*n_actions*n_alphas*n_actions + a_bank*n_actions*n_alphas + s_alpha*n_actions + a_alpha],
                reward + PENALTY - Q[s_bank*n_actions*n_alphas*n_actions + a_bank*n_actions*n_alphas + s_alpha*n_actions + a_alpha],
                Q[s_bank*n_actions*n_alphas*n_actions + a_bank*n_actions*n_alphas + s_alpha*n_actions + a_alpha]/(reward + PENALTY - Q[s_bank*n_actions*n_alphas*n_actions + a_bank*n_actions*n_alphas + s_alpha*n_actions + a_alpha]));
                */

                // UPDATE Q MATRIX
                if (it == max_steps){
                    Q[s_bank*n_actions*n_alphas*n_actions + a_bank*n_actions*n_alphas + s_alpha*n_actions + a_alpha] +=
                        Alpha*(reward - Q[s_bank*n_actions*n_alphas*n_actions + a_bank*n_actions*n_alphas + s_alpha*n_actions + a_alpha]);
                }
                else {
                    Q[s_bank*n_actions*n_alphas*n_actions + a_bank*n_actions*n_alphas + s_alpha*n_actions + a_alpha] +=
                        + Alpha*(reward + Gamma*Q[s_bank1*n_actions*n_actions*n_alphas + a_bank1*n_actions*n_alphas + s_alpha1*n_actions + a_alpha1] \
                        - Q[s_bank*n_actions*n_alphas*n_actions + a_bank*n_actions*n_alphas + s_alpha*n_actions + a_alpha]);
                }

                Q_count[s_bank*n_actions*n_alphas*n_actions + a_bank*n_actions*n_alphas + s_alpha*n_actions + a_alpha] += 1;

                // MOVE ON: S = S1, A = A1

                s_alpha = s_alpha1;
                s_bank = s_bank1;

                a_alpha = a_alpha1;
                a_bank = a_bank1;
                
            }

            it += 1;
        }

        episode += 1;

    }

    //print_mat(Q);

    //printf("save matrix\n");
    //save_matrix(Q, "hcQ_matrix3d.dat");

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
    free(Q_count);

    fclose(out);
    fclose(rew);
    fclose(Q_mat);
    fclose(Q_mat_count);
    fclose(policy);

    return 0;

    }

#include "../../../Dynamics/dynamics_2d_cartesian.h"
#include "../../../Dynamics/winds.h"
#include "sarsa_alpha_vrel_angle.h"

double reward_dt = h*decision_time;
#define dim 2

int main(int argc, char *argv[]){

    if (num_saved_matrices > learning_episodes){
        printf("===> Error: check save matrix step\n");
        return 1;
    }

    int save_matrix_step = (int)learning_episodes/num_saved_matrices;

    FILE *out ,*rew, *Q_mat, *Q_mat_count, *infos;
    out = fopen("cout_epsilon02_8.txt", "w");
    rew = fopen("crewards_epsilon02_8.txt", "w");
    Q_mat = fopen("cQ_mat_epsilon02_8.txt", "w");
    Q_mat_count = fopen("cQ_count_epsilon02_8.txt", "w");
    infos = fopen("infos_epsilon02_8.txt", "w");

    fprintf(rew, "episode,epsilon,learing_rate,steps,return,with_penalty\n");
    fprintf(out, "episode,step,x_kite,z_kite,x_block,z_block,theta,vtheta,windx,windy,v_block,vrelx,vrely,v_rel_mod,alpha,action,vrel_angle_idx,reward,Q[s+0],Q[s+1],Q[s+2]\n");
    fprintf(Q_mat, "episode,vrel_angle_idx,alpha_idx,action_0,action_1,action_2\n");
    fprintf(Q_mat_count, "episode,vrel_angle_idx,alpha_idx,action_0,action_1,action_2\n");
    fprintf(infos, "relative velocities partition\n");

    for (int i=0; i<n_angles+1; i++){
        fprintf(infos,"%f, ", angles[i]);
    }
    fprintf(infos,"\n");
    fclose(infos);

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

    double W[dim] = {0, 0};
    int et_val = 0;

    //double lr = 1.0; // adaptive learning rate
    double vrel_x, vrel_z;
    double vrel_angle;

    double T = 0;
    double F_attr = 0;
    double lift=0, drag=0;
    int sector = 0;
    double l0, l1;
    double d0, d1;

    // ======== LEARNING VARIABLES =======

    double * Q = (double*) malloc(n_alphas * n_actions * n_angles * sizeof(double));
    int * Q_count = (int*) malloc(n_alphas * n_actions * n_angles * sizeof(int));

    double * x_block_seq = (double*) malloc(max_steps * sizeof(double));

    double reward = 0.;
    double tot_reward = 0.;

    // SARSA STATE
    int s_alpha, s_alpha1;
    int s_vrel_angle, s_vrel_angle1;

    // SARSA ACTION
    int a_alpha = 0, a_alpha1 = 0;           

    int it = 0;
    int episode = 0;

    rk[0] = 0;
    rk[1] = R;

    streamfunction2d_hard(rk, W);
    printf("W0 = %f, W1 = %f\n", W[0], W[1]);
    double max_return = sqrt(W[0]*W[0] + W[1]*W[1])*max_steps*h;
    printf("reward max %f\n\n", max_return);

    initialize_Q(Q, 2500);

    //read old Q matrx for each velocity index
    //load_matrix(Q, "last_mat_epsilon02_2.txt");
    print_matrix(Q);

    initialize_Q_count(Q_count, 0);

    printf("Total number of episodes: %d\n", learning_episodes);

    // ======================= STARTING EPISODES ==========================

    while (episode < learning_episodes){

         if (episode == 300000){
            printf("Decreasing learning rate: %f\n", Alpha);
            epsilon = 0.01;
        }
        if (episode == 400000){
            printf("Decreasing learning rate: %f\n", Alpha);
            Alpha = 0.1;
        }
        if (episode == 500000){
            epsilon = 0.001;
            printf("Last episode: greedy\n");
        }
        if (episode == 100000){
            epsilon = 0.0001;
            printf("Last episode: greedy\n");
        }
        if (episode == 120000){
            Alpha = 0.01;
            printf("Last episode: greedy\n");
        }
        if (episode == (int)(learning_episodes - 150)){
            epsilon = 0.0;
            printf("Last episode: greedy\n");
        }
    
        // ======================= EPISODE INITIALIZATION ==========================

        printf("Cartesian code\n");
        printf("Episode: %d\nepsilon: %f\n", episode, epsilon);
        printf("Learning rate: %f\n", Alpha);

        theta = theta0;
        vtheta = vtheta0;
	
        variables_initialization(rk, vk, ak, theta, vtheta, r_block, v_block, a_block, r_diff, v_diff, a_diff);

        streamfunction2d_hard(rk, W);

        s_alpha = s_alpha0;

        vrel_x = vk[0] - W[0];
        vrel_z = vk[1] - W[1];

        vrel_angle = atan2(vrel_z, vrel_x);

        s_vrel_angle = find_state_angle( vrel_angle );

        printf("vrelx=%f, vrelz=%f, angle=%f, angle state=%d\n", vrel_x, vrel_z, vrel_angle, s_vrel_angle );
 
        a_alpha = select_alpha_action(epsilon, Q, s_alpha, s_vrel_angle);

        it = 0;
        reward = 0;
        tot_reward = 0.; 

        printf("theta0 = %f, vel_theta0 = %f\n", theta, vtheta);
        printf("initial s_alpha=%d\n", s_alpha);
        printf("W[0]=%f, W[1]=%f\n", W[0], W[1]);

        integration_trajectory(rk, vk, ak, r_block, v_block, a_block, r_diff, v_diff, a_diff, \
                            &theta, &vtheta, s_alpha, W, &lift, &drag, &T, &F_attr, it, &sector, &l0, &l1, &d0, &d1, &et_val);

        streamfunction2d_hard(rk, W);
        
        while (rk[1] > 0){

            //x_block_seq[it] = r_block[0];

            if ((episode%save_matrix_step == 0 || episode == learning_episodes-1) && it%decision_time == 0){

                fprintf(out, "%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%d,%d,%f,%f,%f,%f\n", episode, it, rk[0], rk[1], 
                    r_block[0], r_block[1], theta, vtheta, W[0], W[1], v_block[0], vrel_x, vrel_z, \
                    sqrt(vrel_x*vrel_x + vrel_z*vrel_z),alphas[s_alpha], a_alpha, s_vrel_angle, reward, \
                    Q[s_vrel_angle*n_alphas*n_actions + s_alpha*n_actions + 0], \
                    Q[s_vrel_angle*n_alphas*n_actions + s_alpha*n_actions + 1], \
                    Q[s_vrel_angle*n_alphas*n_actions + s_alpha*n_actions + 2]);
            }  

            if (it > max_steps){
                printf("MAX STEPS, %d, exiting\n", max_steps);
                printf("final alpha=%d\n", s_alpha);
                printf("return=%f, space percurred=%f\n\n", tot_reward, r_block[0]);

                fprintf(rew, "%d,%f,%f,%d,%f,%f\n", episode, epsilon, Alpha, it, tot_reward, 0.);

                if ( episode%save_matrix_step == 0 || episode == learning_episodes-1 ){
                    fill_Q_mat(Q_mat, Q, episode);
                    fill_Q_count(Q_mat_count, Q_count, episode);
                }

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
            
            // SOME CHECKS
            /*if (check(s_alpha, a_alpha, s_vrel_angle) == 1) {
                printf("CHECK: s_alpha %d\n", s_alpha);
                printf("Matrix index error!!\n");
                episode = learning_episodes;
                break;
            }

            if (a_alpha > 2 || a_alpha < 0){
                printf("a_alpha:%d\n", a_alpha1);
                printf("ERROR IN UPDATE STATE");
                episode = learning_episodes;
                break;
            }*/

            if (rk[1] <= 0.) {
                
                fprintf(rew, "%d,%f,%f,%d,%f,%f\n", episode, epsilon, Alpha, it, tot_reward, PENALTY);

                printf("Kite fallen: z<0, steps=%d, break\n", it);
                printf("final alpha=%d\n", s_alpha);
                printf("return-penalty=%f, space percurred=%f\n\n", tot_reward + PENALTY, r_block[0]);

                if ((episode%save_matrix_step == 0 || episode == learning_episodes-1)){
                    fprintf(out, "%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%d,%d,%f,%f,%f,%f\n", episode, it, rk[0], rk[1], 
                    r_block[0], r_block[1], theta, vtheta, W[0], W[1], v_block[0], vrel_x, vrel_z, \
                    sqrt(vrel_x*vrel_x + vrel_z*vrel_z),alphas[s_alpha], a_alpha, s_vrel_angle, reward, \
                    Q[s_vrel_angle*n_alphas*n_actions + s_alpha*n_actions + 0], \
                    Q[s_vrel_angle*n_alphas*n_actions + s_alpha*n_actions + 1], \
                    Q[s_vrel_angle*n_alphas*n_actions + s_alpha*n_actions + 2]);
                }

                /*if (Q_count[s_vrel_angle*n_alphas*n_actions + s_alpha*n_actions + a_alpha] > 0){
                    lr = 1.0/sqrt(Q_count[s_vrel_angle*n_alphas*n_actions + s_alpha*n_actions + a_alpha]);
                }
                else { lr = 1; }*/

                Q[s_vrel_angle*n_alphas*n_actions + s_alpha*n_actions + a_alpha] += Alpha*(reward + PENALTY  
                                    - Q[s_vrel_angle*n_alphas*n_actions + s_alpha*n_actions + a_alpha]);

                Q_count[s_vrel_angle*n_alphas*n_actions + s_alpha*n_actions + a_alpha] += 1;

                if ( episode%save_matrix_step == 0 || episode == learning_episodes-1 ){
                    fill_Q_mat(Q_mat, Q, episode);
                    fill_Q_count(Q_mat_count, Q_count, episode);
                }
                break;
            }

            if (it%decision_time == 0){

                /*if (Q_count[s_vrel_angle*n_alphas*n_actions + s_alpha*n_actions + a_alpha] > 0){
                    lr = 1.0/sqrt(Q_count[s_vrel_angle*n_alphas*n_actions + s_alpha*n_actions + a_alpha]);
                }
                else { lr = 1; }*/

                // UPDATE Q MATRIX

                if (it == max_steps){

                    Q[s_vrel_angle*n_alphas*n_actions + s_alpha*n_actions + a_alpha] += Alpha*(reward
                        - Q[s_vrel_angle*n_alphas*n_actions + s_alpha*n_actions + a_alpha]);
                }
                else {

                    // FIND STATE S1 USING ACTION A

                    s_alpha1 = update_state(s_alpha, a_alpha);

                    vrel_x = vk[0] - W[0];
                    vrel_z = vk[1] - W[1];
                    vrel_angle = atan2(vrel_z, vrel_x);
                    s_vrel_angle1 = find_state_angle(vrel_angle);
                    
                    // SEARCH FOR NEXT ACTION A1

                    a_alpha1 = select_alpha_action(epsilon, Q, s_alpha1, s_vrel_angle1);

                    if (episode%100 == 0){
                        printf("%f, %f, %f, %f\n", Q[s_vrel_angle*n_alphas*n_actions + s_alpha*n_actions + a_alpha],
                        reward + Gamma*Q[s_vrel_angle1*n_alphas*n_actions + s_alpha1*n_actions + a_alpha1]
                        - Q[s_vrel_angle*n_alphas*n_actions + s_alpha*n_actions + a_alpha],
                        Q[s_vrel_angle*n_alphas*n_actions + s_alpha*n_actions + a_alpha]/(reward  
                        + Gamma*Q[s_vrel_angle1*n_alphas*n_actions + s_alpha1*n_actions + a_alpha1]
                        - Q[s_vrel_angle*n_alphas*n_actions + s_alpha*n_actions + a_alpha]), \
                        (reward + Gamma*Q[s_vrel_angle1*n_alphas*n_actions + s_alpha1*n_actions + a_alpha1]
                        - Q[s_vrel_angle*n_alphas*n_actions + s_alpha*n_actions + a_alpha])/Q[s_vrel_angle*n_alphas*n_actions + s_alpha*n_actions + a_alpha]);
                    }

                    Q[s_vrel_angle*n_alphas*n_actions + s_alpha*n_actions + a_alpha] += Alpha*(reward  
                        + Gamma*Q[s_vrel_angle1*n_alphas*n_actions + s_alpha1*n_actions + a_alpha1]
                        - Q[s_vrel_angle*n_alphas*n_actions + s_alpha*n_actions + a_alpha]);  

                    // MOVE ON: S = S1, A = A1
                    s_alpha = s_alpha1;
                    s_vrel_angle = s_vrel_angle1;

                    a_alpha = a_alpha1;
                }

                Q_count[s_vrel_angle*n_alphas*n_actions + s_alpha*n_actions + a_alpha] += 1;
                
            }

            it += 1;

        }

        episode += 1;

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
    free(x_block_seq);
    free(Q_count);

    fclose(out);
    fclose(rew);
    fclose(Q_mat);
    fclose(Q_mat_count);
    fclose(infos);

    return 0;

    }

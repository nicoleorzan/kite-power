#include "../../../Dynamics/dynamics_2d_spherical_cramer.h"
#include "../../../Dynamics/winds.h"
#include "sarsa_alpha_theta_vtheta.h"

#define decision_time 1000

double reward_dt = h*decision_time;

int main(int argc, char *argv[]){

    if (save_matrix_step > learning_episodes){
        printf("===> Error: check save matrix step\n");
        return 1;
    }

    FILE *out ,*rew, *Q_mat, *Q_mat_count, *policy;
    out = fopen("pout.txt", "w");
    rew = fopen("prewards.txt", "w");
    Q_mat = fopen("pQ_mat.txt", "w");
    Q_mat_count = fopen("pQ_count.txt", "w");
    policy = fopen("ppolicy.txt", "w");
    fprintf(rew, "episode,epsilon,Alpha,steps,return\n");
    fprintf(out, "t         x_kite          z_kite         r_block          z_blocco          wind_x       wind_y       v_blocco_x\n");
    fprintf(Q_mat, "episode,alpha_idx,theta,vtheta,action_0,action_1,action_2\n");
    fprintf(policy, "step        alpha       action      reward        Q[s+0]      Q[s+1]      Q[s+2]\n");
    
    // ======== DYNAMICS VARIABLES =======

    // vettori moto kite dall'origine fissa (x, z)
    double *rk = (double*) malloc(dim * sizeof(double)); 
    double *vk = (double*) malloc(dim * sizeof(double)); 
    double *ak = (double*) malloc(dim * sizeof(double)); 

    // vettori moto blocco dall'origine fissa (x, z)
    double *r_block = (double*) malloc(dim * sizeof(double)); 
    double *v_block = (double*) malloc(dim * sizeof(double));
    double *a_block = (double*) malloc(dim * sizeof(double));  
    
    // theta, dtheta, ddtheta
    double *theta = (double*) malloc(3 * sizeof(double)); 

    double W[dim];

    double T = 0;
    double F_attr = 0;
    double lift=0, drag=0;
    int sector = 0; 

    // ======== LEARNING VARIABLES =======

    double * Q = (double*) malloc(n_alphas * n_actions * n_thetas * n_thetas_vel * sizeof(double));
    int * Q_count = (int*) malloc(n_alphas * n_actions  * n_thetas * n_thetas_vel * sizeof(int));

    double reward = 0;
    double tot_reward = 0;

    // SARSA STATE
    int s_alpha, s_alpha1;
    int s_theta, s_theta1;
    int s_vtheta, s_vtheta1;

    // SARSA ACTION
    int a_alpha = 0, a_alpha1 = 0;           

    int it = 0;
    int episode = 0;

    rk[0] = 0;
    rk[1] = R;
    streamfunction2d(rk, W);
    printf("W0 = %f, W1 = %f\n", W[0], W[1]);
    double reward_max = sqrt(W[0]*W[0] + W[1]*W[1])*max_steps*h;
    printf("reward max %f\n\n", reward_max);

    initialize_Q(Q, reward_max);
    initialize_Q_count(Q_count, 0);

    printf("Total number of episodes: %d\n", learning_episodes);

    // ======================= STARTING EPISODES ==========================

    while (episode < learning_episodes){

        if (episode == (int)(learning_episodes/4)){
            Alpha = Alpha*0.1;
            printf("Decreasing learning rate: %f\n", Alpha);
        }
        if (episode == (int)(learning_episodes/2)){
            epsilon = epsilon + 0.05;
            Alpha = Alpha*0.1;
            printf("Decreasing learning rate: %f\n", Alpha);
        }
        if (episode == (int)(learning_episodes*3/4)){
            Alpha = Alpha*0.1;
            printf("Decreasing learning rate: %f\n", Alpha);
        }
    
        // ======================= EPISODE INITIALIZATION ==========================

        printf("Episode: %d\nepsilon: %f\n", episode, epsilon);
        printf("Learning rate: %f\n", Alpha);

        theta[0] = theta0; 
        theta[1] = vtheta0;
	
        variables_initialization(rk, vk, ak, theta0, vtheta0, r_block, v_block, a_block, theta);
        
        streamfunction2d(rk, W);

        s_alpha = s_alpha0;     
        s_theta = find_state(theta[0], thetas);
        s_vtheta = find_state(theta[1], thetas_vel);
        a_alpha = select_alpha_action(epsilon, Q, s_alpha, s_theta, s_vtheta, episode);

        it = 0;
        reward = 0;
        tot_reward = 0.; 

        printf("theta0 = %f, vel_theta0 = %f\n", theta[0], theta[1]);
        printf("initial alpha state=%d\n", s_alpha);
        printf("W[0]=%f\n", W[0]);
        printf("W[1]=%f\n", W[1]);

        integration_trajectory(rk, vk, ak, r_block, v_block, a_block, theta, s_alpha, \
                             W, &lift, &drag, &T, &F_attr, it, &sector);

        streamfunction2d(rk, W);
        
        while (rk[1] > 0){

            if (episode == learning_episodes - 1 && it%decision_time == 0){
                fprintf(policy, "%d      %f     %d     %f     %f     %f     %f\n", \
                it, alphas[s_alpha], a_alpha, reward, 
                Q[s_alpha*n_actions*n_thetas*n_thetas_vel + 0*n_thetas*n_thetas_vel 
                + s_theta*n_thetas_vel + s_vtheta],
                Q[s_alpha*n_actions*n_thetas*n_thetas_vel + 1*n_thetas*n_thetas_vel
                 + s_theta*n_thetas_vel + s_vtheta], 
                Q[s_alpha*n_actions*n_thetas*n_thetas_vel + 2*n_thetas*n_thetas_vel
                 + s_theta*n_thetas_vel + s_vtheta]);

                fprintf(out, "%d       %f       %f      %f      %f      %f      %f     %f\n", it, \
                    rk[0], rk[1], r_block[0], r_block[1], W[0], W[1], v_block[0]);
            }

            if (it > max_steps){
                printf("MAX STEPS, %d, exiting\n", max_steps);
                printf("return=%f, space percurred=%f\n\n", tot_reward, r_block[0]);
                fprintf(rew, "%d,%f,%f,%d,%f\n", episode, epsilon, Alpha, it, tot_reward);

                if ( episode%save_matrix_step == 0){
                    fill_Q_mat(Q_mat, Q, episode);
                    fill_Q_count(Q_mat_count, Q_count);
                }

                break;
            }

            integration_trajectory(rk, vk, ak, r_block, v_block, a_block, theta, s_alpha, \
                             W, &lift, &drag, &T, &F_attr, it, &sector);
            
            streamfunction2d(rk, W);

            reward = fabs(v_block[0])*reward_dt;

            if (it%decision_time == 0){
                tot_reward += reward;
            }
            
            // SOME CHECKS
            /*if (check(s_alpha, a_alpha, s_theta) == 1) {
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
                
                fprintf(rew, "%d,%f,%f,%d,%f\n", episode, epsilon, Alpha, it, tot_reward);

                printf("Kite fallen: z<0, steps=%d, break\n", it);
                printf("return=%f, space percurred=%f\n\n", tot_reward, r_block[0]);

                Q[s_alpha*n_actions*n_thetas*n_thetas_vel + a_alpha*n_thetas*n_thetas_vel 
                    + s_theta*n_thetas_vel + s_vtheta] += Alpha*(reward + PENALTY 
                    - Q[s_alpha*n_actions*n_thetas*n_thetas_vel + a_alpha*n_thetas*n_thetas_vel 
                    + s_theta*n_thetas_vel + s_vtheta]);

                Q_count[s_alpha*n_actions*n_thetas*n_thetas_vel + a_alpha*n_thetas*n_thetas_vel 
                    + s_theta*n_thetas_vel + s_vtheta] += 1;

                if ( episode%save_matrix_step == 0){
                    fill_Q_mat(Q_mat, Q, episode);
                    fill_Q_count(Q_mat_count, Q_count);
                }
                break;
            }

            if (it%decision_time == 0){

                // FIND STATE S1 USING ACTION A

                s_alpha1 = update_state(s_alpha, a_alpha);
                s_theta1 = find_state(theta[0], thetas);
                s_vtheta1 = find_state(theta[1], thetas_vel);

                // SEARCH FOR NEXT ACTION A1

                a_alpha1 = select_alpha_action(epsilon, Q, s_alpha1, s_theta1, s_vtheta1, episode);

                // UPDATE Q MATRIX
                if (it == max_steps){

                    Q[s_alpha*n_actions*n_thetas*n_thetas_vel + a_alpha*n_thetas*n_thetas_vel 
                        + s_theta*n_thetas_vel + s_vtheta] += Alpha*(reward
                        - Q[s_alpha*n_actions*n_thetas*n_thetas_vel + a_alpha*n_thetas*n_thetas_vel 
                        + s_theta*n_thetas_vel + s_vtheta]);
                }
                else {

                    Q[s_alpha*n_actions*n_thetas*n_thetas_vel + a_alpha*n_thetas*n_thetas_vel 
                        + s_theta*n_thetas_vel + s_vtheta] += Alpha*(reward
                        + Gamma*Q[s_alpha1*n_actions*n_thetas*n_thetas_vel + a_alpha1*n_thetas*n_thetas_vel 
                        + s_theta1*n_thetas_vel + s_vtheta1]
                        - Q[s_alpha*n_actions*n_thetas*n_thetas_vel + a_alpha*n_thetas*n_thetas_vel 
                        + s_theta*n_thetas_vel + s_vtheta]);  
                }

                Q_count[s_alpha*n_actions*n_thetas*n_thetas_vel + a_alpha*n_thetas*n_thetas_vel 
                    + s_theta*n_thetas_vel + s_vtheta] += 1;

                // MOVE ON: S = S1, A = A1

                s_alpha = s_alpha1;
                s_theta = s_theta1;
                s_vtheta = s_vtheta1;

                a_alpha = a_alpha1;
                
            }

            it += 1;
        }

        episode += 1;

    }

    //print_mat(Q);

    //printf("save matrix\n");
    //save_matrix(Q, "Q_matrix.dat");
    //printf("load matrix\n");
    //load_matrix(Q2, "Q_matrix.dat");
    //print_mat(Q2);

    free(rk);
    free(vk);
    free(ak);

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

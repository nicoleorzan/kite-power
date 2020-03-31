#include "../../../Dynamics/dynamics_2d_spherical_cramer.h"
#include "../../../Dynamics/winds.h"
#include "sarsa_alpha.h"

double reward_dt = h*decision_time;

int main(int argc, char *argv[]){

    if (num_saved_matrices > learning_episodes){
        printf("===> Error: check save matrix step\n");
        return 1;
    }

    int save_matrix_step = (int)learning_episodes/num_saved_matrices;

    FILE *out ,*rew, *Q_mat, *Q_mat_count, *policy;
    out = fopen("pout.txt", "w");
    rew = fopen("prewards.txt", "w");
    Q_mat = fopen("pQ_matrix.txt", "w");
    Q_mat_count = fopen("pQ_counter.txt", "w");
    policy = fopen("ppolicy.txt", "w");
    fprintf(rew, "episode,epsilon,Alpha,steps,return\n");
    fprintf(out, "t,x_kite,z_kite,x_block,z_block,theta,vtheta,windx,windy,v_block,Tension\n");
    fprintf(Q_mat, "episode,alpha_idx,action_0,action_1,action_2\n");
    fprintf(Q_mat_count, "episode,alpha_idx,action_0,action_1,action_2\n");
    fprintf(policy, "step,alpha,action,reward,Q[s+0],Q[s+1],Q[s+2]\n");

    // ======== DYNAMICS VARIABLES =======

    // vettori moto kite dall'origine fissa (x, z)
    double *rk = (double*) malloc(dim * sizeof(double)); 
    double *vk = (double*) malloc(dim * sizeof(double)); 
    double *ak = (double*) malloc(dim * sizeof(double)); 

    // vettori moto block dall'origine fissa (x, z)
    double *r_block = (double*) malloc(dim * sizeof(double)); 
    double *v_block = (double*) malloc(dim * sizeof(double));
    double *a_block = (double*) malloc(dim * sizeof(double));

    // theta, dtheta, ddtheta
    double *theta = (double*) malloc(3 * sizeof(double));  

    double T = 0;
    double F_attr = 0;
    double lift = 0;
    double drag = 0;
    int sector = 0;

    // ======== LEARNING VARIABLES =======

    double * Q = (double*) malloc(n_alphas * n_actions * sizeof(double));
    int * Q_count = (int*) malloc(n_alphas * n_actions * sizeof(int));

    double reward = 0;
    double tot_reward = 0;

    // SARSA STATE
    int s_alpha, s_alpha1;

    // SARSA ACTION
    int a_alpha = 0, a_alpha1 = 0;           

    int it = 0;
    int episode = 0;

    double W[dim] = {20, 0};

    //rk[0] = 0;
    //rk[1] = R;

    //streamfunction2d(rk, W);
    printf("W0 = %f, W1 = %f\n", W[0], W[1]);
    double reward_max = sqrt(W[0]*W[0] + W[1]*W[1])*max_steps*h;
    printf("reward max %f\n\n", reward_max);

    initialize_Q(Q, reward_max);
    initialize_Q_count(Q_count, 0);

    //load_matrix(Q, "Q_matrix.dat");

    printf("Total number of episodes: %d\n", learning_episodes);

    // ======================= STARTING EPISODES ==========================

    while (episode < learning_episodes){

<<<<<<< HEAD
<<<<<<< HEAD
        if (episode == (int)(learning_episodes/2)){
=======
        /*if (episode == (int)(learning_episodes/5)){
>>>>>>> 8a4af4f473d680572248ff2f1a7b5406024ab1e4
            Alpha = Alpha*0.1;
            printf("Decreasing learning rate: %f\n", Alpha);
        }
        if (episode == (int)(learning_episodes*3/4)){
            epsilon += 0.05;
=======
        /*if (episode == (int)(learning_episodes/5)){
            Alpha = Alpha*0.1;
            printf("Decreasing learning rate: %f\n", Alpha);
<<<<<<< HEAD
        }
        if (episode == (int)(learning_episodes*2/5)){
            Alpha = Alpha*0.1;
            printf("Decreasing learning rate: %f\n", Alpha);
=======
>>>>>>> 8a4af4f473d680572248ff2f1a7b5406024ab1e4
        }*/
        if (episode == (int)(learning_episodes*3/5)){
            epsilon = epsilon + 0.05;
            Alpha = Alpha*0.1;
            printf("Decreasing learning rate: %f\n", Alpha);
        }
        if (episode == (int)(learning_episodes*5.4)){
>>>>>>> 8a4af4f473d680572248ff2f1a7b5406024ab1e4
            Alpha = Alpha*0.1;
            printf("Decreasing learning rate: %f\n", Alpha);
        }
    
        // ======================= EPISODE INITIALIZATION ==========================

        printf("Polar code\n");
        printf("Episode: %d\nepsilon: %f\n", episode, epsilon);
        printf("Learning rate: %f\n", Alpha);

        theta[0] = theta0;
        theta[1] = vtheta0;

        variables_initialization(rk, vk, ak, theta[0], theta[1], r_block, v_block, a_block, theta);

        //streamfunction2d(rk, W);

        s_alpha = s_alpha0;
        a_alpha = select_alpha_action(epsilon, Q, s_alpha, episode);

        it = 0;
        reward = 0;
        tot_reward = 0.; 

        printf("theta0 = %f, vel_theta0 = %f\n", theta[0], theta[1]);
        printf("initial alpha state=%d\n", s_alpha);
        printf("W[0]=%f\n", W[0]);
        printf("W[1]=%f\n", W[1]);

        integration_trajectory(rk, vk, ak, r_block, v_block, a_block, theta, s_alpha, \
                             W, &lift, &drag, &T, &F_attr, it, &sector);

        //streamfunction2d(rk, W);
        
        while (rk[1] > 0){

            if (episode == learning_episodes - 1 && it%decision_time == 0){
                fprintf(policy, "%d,%f,%d,%f,%f,%f,%f\n", \
                it, alphas[s_alpha], a_alpha, reward, Q[s_alpha*n_actions + 0], \
                Q[s_alpha*n_actions + 1], Q[s_alpha*n_actions + 2]);
            }

            if (it > max_steps){
                printf("MAX STEPS, %d, exiting\n", max_steps);
                printf("return=%f, space percurred=%f\n\n", tot_reward, r_block[0]);

                fprintf(rew, "%d,%f,%f,%d,%f\n", episode, epsilon, Alpha, it, tot_reward);

                if ( episode%save_matrix_step == 0){
                    fill_Q_mat(Q_mat, Q, episode);
                }
                if (episode == learning_episodes-1){
                    fill_Q_count(Q_mat_count, Q_count, episode);
                }

                break;
            }

            integration_trajectory(rk, vk, ak, r_block, v_block, a_block, theta, s_alpha, \
                             W, &lift, &drag, &T, &F_attr, it, &sector);

            //streamfunction2d(rk, W);

            reward = fabs(v_block[0])*reward_dt;

            if (it%decision_time == 0){
                tot_reward += reward;
            }
            
            // control if alpha index is bigger or equal then 15 or smaller than 0
            /*if (check(s_alpha) == 1) {
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
            
            // save data of last episode for plot

            if ( (episode == learning_episodes - 1) && (it%decision_time == 0) ){
                fprintf(out, "%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", it, rk[0], rk[1], 
                    r_block[0], r_block[1], theta[0], theta[1], W[0], W[1], v_block[0], T);
            }

            if (rk[1] <= 0.) {

                fprintf(rew, "%d,%f,%f,%d,%f\n", episode, epsilon, Alpha, it, tot_reward);

                printf("Kite fallen: z<0, steps=%d, break\n", it);
                printf("return=%f, space percurred=%f\n\n", tot_reward, r_block[0]);

                Q[s_alpha*n_actions + a_alpha] = Q[s_alpha*n_actions + a_alpha] + \
                        Alpha*(reward + PENALTY - Q[s_alpha*n_actions + a_alpha]);

                Q_count[s_alpha*n_actions + a_alpha] += 1;

                if ( episode%save_matrix_step == 0){
                    fill_Q_mat(Q_mat, Q, episode);
                }
                if (episode == learning_episodes-1){
                    fill_Q_count(Q_mat_count, Q_count, episode);
                }
<<<<<<< HEAD
<<<<<<< HEAD
                
=======
=======
>>>>>>> 8a4af4f473d680572248ff2f1a7b5406024ab1e4

>>>>>>> 8a4af4f473d680572248ff2f1a7b5406024ab1e4
                break;
            }

            if (it%decision_time == 0){

                // FIND STATE S1 USING ACTION A

                s_alpha1 = update_state(s_alpha, a_alpha);

                // SEARCH FOR NEXT ACTION A1

                a_alpha1 = select_alpha_action(epsilon, Q, s_alpha1, episode);

                // UPDATE Q MATRIX
                if (it == max_steps){
                    Q[s_alpha*n_actions + a_alpha] = Q[s_alpha*n_actions + a_alpha] \
                            + Alpha*(reward - Q[s_alpha*n_actions + a_alpha]);
                }
                else {
                    Q[s_alpha*n_actions + a_alpha] = Q[s_alpha*n_actions + a_alpha] \
                            + Alpha*(reward + Gamma*Q[s_alpha1*n_actions + a_alpha1] \
                            - Q[s_alpha*n_actions + a_alpha]);
                }

                Q_count[s_alpha*n_actions + a_alpha] += 1;

                // MOVE ON: S = S1, A = A1

                s_alpha = s_alpha1;
                a_alpha = a_alpha1;
                
            }

            it += 1;
        }

        episode += 1;

    }

    print_mat(Q);

    printf("save matrix\n");
    save_matrix(Q, "pQ_matrix.dat");
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

    free(theta);

    fclose(out);
    fclose(rew);
    fclose(Q_mat);
    fclose(Q_mat_count);
    fclose(policy);

    return 0;

    }

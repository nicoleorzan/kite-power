#include "../../Dynamics/dynamics_2d.h"
#include "sarsa_alpha.h"

#define decision_time 1000

#define theta0 0.
#define vtheta0 0.1
#define k 0.5

double reward_dt = h*decision_time;

int main(int argc, char *argv[]){

    FILE *out ,*rew, *Q_mat, *policy;
    out = fopen("out1", "w");
    rew = fopen("rewards1", "w");
    Q_mat = fopen("Q_matrix1", "w");
    policy = fopen("policy1", "w");
    fprintf(out, "t         x_kite          z_kite         x_blocco          z_blocco          wind_x       wind_y       v_blocco_x\n");
    fprintf(out, "t       Q1      Q2      Q3      Q4      Q5      Q6      Q7      Q8      Q9      Q10     Q11     Q12      Q13     Q14n");
    fprintf(policy, "step        alpha       action      reward        Q[s+0]      Q[s+1]      Q[s+2]\n");

    // vettori moto kite dall'origine fissa (x, z)
    double *rk = (double*) malloc(2 * sizeof(double)); 
    double *vk = (double*) malloc(2 * sizeof(double)); 
    double *ak = (double*) malloc(2 * sizeof(double)); 

    // vettori moto blocco dall'origine fissa (x, z)
    double *x_blocco = (double*) malloc(2 * sizeof(double)); 
    double *v_blocco = (double*) malloc(2 * sizeof(double));
    double *a_blocco = (double*) malloc(2 * sizeof(double));  
    
    // theta, dtheta, ddtheta
    double *theta = (double*) malloc(3 * sizeof(double));  
    double power = 0;
    double T = 0;

    double W_max = 50*k;

    double * Q = (double*) malloc(n_alphas * n_actions * sizeof(double));

    double reward = 0;
    double tot_reward = 0;
    double penalty = -200;
    double Cl, Cd;

    double theta_star;

    double Wx, Wy;
    Wx = 30;
    Wy = 0;
    double lift=0, drag=0;
    double fluct = 0;

    // SARSA STATE
    int s_alpha, s_alpha1;

    // SARSA ACTION
    int a_alpha = 0, a_alpha1 = 0;           

    int episode = 0;
    double epsilon;
    int vv = 0;

    double reward_max = Wx*max_steps*h;
    
    printf("W max %f\n\n", W_max);
    printf("reward max %f\n\n", reward_max);

    initialize_Q(Q, reward_max);

    //load_matrix(Q, "Q_matrix.dat");
    print_mat(Q);
    //print_mat(Q);

    printf("Total number of episodes: %d\n", learning_episodes);

    // ======================= STARTING EPISODES ==========================

    while (episode < learning_episodes){

        //epsilon = (double)episode/(double)learning_episodes;
        epsilon = 0.9;

        if (episode == (int)(learning_episodes/2)){
            Alpha = Alpha*0.1;
            printf("Decreasing learning rate: %f\n", Alpha);
        }
        if (episode == (int)(learning_episodes*3/4)){
            Alpha = Alpha*0.1;
            printf("Decreasing learning rate: %f\n", Alpha);
        }
    
        // ======================= EPISODE INITIALIZATION ==========================

        printf("Episode: %d\nepsilon: %f\n", episode, epsilon);

        variables_initialization(x_blocco, v_blocco, a_blocco, theta, rk, vk, ak, theta0, vtheta0);
        Wx = 20;
        Wy = 10;

        printf("theta0 = %f, vel_theta0 = %f\n", theta[0], theta[1]);

        s_alpha = 10;
        
        Cl = CL_alpha[s_alpha];
        Cd = CD_alpha[s_alpha];

        printf("init alpha index=%d, Cl=%f, Cd=%f\n", s_alpha, Cl, Cd);

        a_alpha = select_alpha_action(epsilon, Q, s_alpha, episode);

        int it = 0;

        reward = 0;
        tot_reward = 0.; 

        //W = 0.6*rk[1];   //rand() % (W_max - W_min) + W_min; // random num between 7 and 40
        printf("Wx=%f\n", Wx);
        printf("Wy=%f\n", Wy);

        integration_trajectory(rk, vk, ak, x_blocco, v_blocco, a_blocco, theta, &T, \
                                   &power, Cl, Cd, Wx, Wy, &lift, &drag);

        if (episode == learning_episodes - 1 ){
            fprintf(policy, "%d      %f     %d     %f     %f     %f     %f\n", \
            it, alphas[s_alpha], a_alpha, reward, Q[s_alpha*n_actions + 0], \
            Q[s_alpha*n_actions + 1], Q[s_alpha*n_actions + 2]);
        }
        
        while (rk[1] > 0){

            if (episode == learning_episodes - 1 && it%decision_time == 0){
                fprintf(policy, "%d      %f     %d     %f     %f     %f     %f\n", \
                it, alphas[s_alpha], a_alpha, reward, Q[s_alpha*n_actions + 0], \
                Q[s_alpha*n_actions + 1], Q[s_alpha*n_actions + 2]);

                /*printf("step        alpha       action      reward        Q[s]      Q[s+1]      Q[s+2]\n");

                printf("%d      %f     %d     %f     %f     %f     %f\n", \
                it, alphas[s_alpha], a_alpha, reward, Q[s_alpha*n_actions + 0], \
                Q[s_alpha*n_actions + 1], Q[s_alpha*n_actions + 2]);*/
            }   

            if (it > max_steps){
                printf("MAX STEPS, %d, exiting\n", max_steps);
                printf("return=%f, space percurred=%f\n\n", tot_reward, x_blocco[0]);

                fprintf(rew, "%d    %f\n", episode, tot_reward);

                fprintf(Q_mat,"%d ", episode);
                for (int i=0; i<n_alphas; i++){
                    for (int j=0; j<n_actions; j++){
                        fprintf(Q_mat,"%f ", Q[i*n_actions + j]);
                    }
                    fprintf(Q_mat,"\n");
                }
                fprintf(Q_mat,"\n");

                tot_reward = 0;
                break;
            }

            integration_trajectory(rk, vk, ak, x_blocco, v_blocco, a_blocco, theta, &T, \
                                   &power, Cl, Cd, Wx, Wy, &lift, &drag);
            //W = wind_const*rk[1];

            theta_star = atan((lift - m*g)/drag); 

            reward = fabs(v_blocco[0])*reward_dt;

            if (it%decision_time == 0){
                tot_reward += reward;
            }
            
            // control if alpha index is bigger or equal then 15 or smaller than 0
            if (check(s_alpha) == 1) {
                printf("CHECK: s_alpha %d\n", s_alpha);
                printf("Matrix index error!!\n");
                episode = learning_episodes;
                break;
            }

            // save data of last episode for plot
            if ( (episode == learning_episodes - 1) && (it%decision_time == 0) ){
                fprintf(out, "%d       %f       %f      %f      %f      %f      %f     %f\n", it, \
                    rk[0], rk[1], x_blocco[0], x_blocco[1], Wx, Wy, v_blocco[0]);
            }

            if (rk[1] <= 0.) {
                fprintf(rew, "%d    %f\n", episode, tot_reward);

                fprintf(Q_mat,"%d ", episode);
                for (int i=0; i<n_alphas; i++){
                    for (int j=0; j<n_actions; j++){
                        fprintf(Q_mat,"%f ", Q[i*n_actions + j]);
                    }
                }
                fprintf(Q_mat,"\n");
                printf("Kite fallen: z<0, steps=%d, break\n", it);
                printf("return=%f, space percurred=%f\n\n", tot_reward, x_blocco[0]);

                Q[s_alpha*n_actions + a_alpha] = Q[s_alpha*n_actions + a_alpha] + \
                        Alpha*(reward + penalty - Q[s_alpha*n_actions + a_alpha]);
                
                tot_reward = 0;
                break;
            }

            if (a_alpha > 2 || a_alpha < 0){
                printf("a_alpha:%d\n", a_alpha1);
                printf("ERROR IN UPDATE STATE");
                episode = learning_episodes;
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

                // MOVE ON: S = S1, A = A1

                s_alpha = s_alpha1;
                a_alpha = a_alpha1;
                
            }

            Cl = CL_alpha[s_alpha];
            Cd = CD_alpha[s_alpha];

            it += 1;
        }

        episode += 1;

    }

    print_mat(Q);

    printf("save matrix\n");
    save_matrix(Q, "Q_matrix.dat");
    printf("load matrix\n");
    //load_matrix(Q2, "Q_matrix.dat");
    //print_mat(Q2);


    free(rk);
    free(vk);
    free(ak);
    free(Q);
    //free(Q2);
    free(x_blocco);
    free(v_blocco);
    free(a_blocco);
    free(theta);

    fclose(out);
    fclose(rew);
    fclose(Q_mat);
    fclose(policy);
    
    remove("a.out");

    return 0;

    }

#include "../../Dynamics/dynamics_3d_cartesian.h"
#include <time.h>
#include <string.h>
#include <stdbool.h>
#include "sarsa-bank-attack.h"

#define decision_time 10000

#define theta0 PI/4.
#define phi0 0.
#define dim 3

double reward_dt = h*decision_time;

void fill_Q_file(FILE *Q_mat, double *Q){
  
  for (int i=0; i<n_alphas; i++){
    for (int j=0; j<n_act; j++){
      for (int kk=0; kk<n_bank; kk++){
	      for (int l=0; l<n_act; l++){
	        fprintf(Q_mat,"%f ", Q[i*n_act*n_bank*n_act + j*n_bank*n_act + kk*n_act + l]);
	      }   
      }
    }
  }
  fprintf(Q_mat,"\n");
  
};

void fill_Q_counter(FILE *Q_mat_counter, int* Q_counter){
  
    //fprintf(Q_mat_counter,"%d ", episode);
    for (int i=0; i<n_alphas; i++){
      for (int j=0; j<n_act; j++){
        for (int kk=0; kk<n_bank; kk++){
          for (int l=0; l<n_act; l++){
          fprintf(Q_mat_counter,"%d ", Q_counter[i*n_act*n_bank*n_act + j*n_bank*n_act + kk*n_act + l]);
          }
	      }
      }
    }
    fprintf(Q_mat_counter,"\n");

};

// ============== FILE INPUT: ATTACK ANGLE AND WIND COEFF AND IF TRAJECTORY IS NEEDED ============

int main(int argc, char *argv[]){

  // ========================= READING INPUT VARIABLES ==========================

  FILE *out ,*rew, *Q_mat, *policy, *Q_mat_counter;
  out = fopen("out", "w");
  rew = fopen("rewards", "w");
  Q_mat = fopen("Q_mat", "w");
  Q_mat_counter = fopen("Q_count", "w");
  policy = fopen("policy", "w");
  
  fprintf(out, "t      x_kite     y_kite     z_kite     x_blocco    y_blocco     z_blocco      wind_x     wind_y    wind_z   v_blocco_x     v_blocco_y\n");
  fprintf(policy, "step        alpha       action      reward        Q[s+0]      Q[s+1]      Q[s+2]\n");
    
  // ============================ VARIABLES DEFINITION ============================

  // kite motion vectors from fixed origin (x, z)
  double *rk = (double*) malloc(dim * sizeof(double)); 
  double *rk1 = (double*) malloc(dim * sizeof(double)); 
  double *vk = (double*) malloc(dim * sizeof(double)); 
  double *ak = (double*) malloc(dim * sizeof(double)); 
  
  // block motion vectors from fixed origin  (x, z)
  double *r_block = (double*) malloc(dim * sizeof(double)); 
  double *v_block = (double*) malloc(dim * sizeof(double));
  double *a_block = (double*) malloc(dim * sizeof(double));  
  
  // block motion vectors from fixed origin  (x, z)
  double *r_diff = (double*) malloc(dim * sizeof(double)); 
  double *v_diff = (double*) malloc(dim * sizeof(double));
  double *a_diff = (double*) malloc(dim * sizeof(double));

  double * Q = (double*) malloc(n_alphas * n_act * n_bank * n_act *  sizeof(double));
  int * Q_counter = (int*) malloc(n_alphas * n_act * n_bank * n_act *  sizeof(int));

  double theta = theta0;
  double phi = phi0;
  double r_diff_modulo;// = sqrt(r_diff[0]*r_diff[0] + r_diff[1]*r_diff[1] + r_diff[2]*r_diff[2]);

  double reward = 0;
  double tot_reward = 0;
  double penalty = -500;
    
  double lift=0, drag=0;
  double F_vinc;
  double theta_star;
  double T = 0;
  int stability = 0;
  int decollato = 0;

  // SARSA STATE
  int s_alpha, s_alpha1;
  int s_mu, s_mu1;

  // SARSA ACTION
  int a_alpha = 0, a_alpha1 = 0;
  int a_mu = 0, a_mu1 = 0;

  int episode = 0, it = 0;

	double W[3] = {10, 0, 0};
	double Wmax = sqrt(W[0]*W[0] + W[1]*W[1] + W[2]*W[2]);
  double reward_max = Wmax*max_steps*h;
  printf("W max %f\n\n", Wmax);
  printf("reward max %f\n\n", reward_max);

  initialize_Q(Q, 4000);
  initialize_Q_counter(Q_counter, 0);

  printf("Total number of episodes: %d\n", learning_episodes);
    
  // ======================= STARTING EPISODES ==========================

  while (episode < learning_episodes){

    if (episode == (int)(learning_episodes/4)){
      Alpha = Alpha*0.1;
      epsilon = 0.85;
      printf("Decreasing learning rate: %f\n", Alpha);
      printf("Decreasing exploration rate (increasing epsilon): %f\n", epsilon);
    }
    if (episode == (int)(learning_episodes/2)){
      Alpha = Alpha*0.1;
      epsilon = 0.90;
      printf("Decreasing learning rate: %f\n", Alpha);
      printf("Decreasing exploration rate (increasing epsilon): %f\n", epsilon);
    }
    if (episode == (int)(learning_episodes*3/4)){
      Alpha = Alpha*0.1;
      epsilon = 0.96;
      printf("Decreasing learning rate: %f\n", Alpha);
      printf("Decreasing exploration rate (increasing epsilon): %f\n", epsilon);
    }
    if (episode == learning_episodes-1){
      epsilon = 0.99;
    }

    // ======================= EPISODE INITIALIZATION ==========================

    printf("EPISODE: %d\nepsilon: %f\n", episode, epsilon);
    printf("Learning rate: %f\n", Alpha);

    s_alpha = 10;
    s_mu = 0;

    printf("initial attack angle=%d, bank_angle=%d\n", s_alpha, s_mu);

    variables_initialization(rk, vk, ak, theta0, phi0, r_block, v_block, a_block);

    printf("theta = %f, phi = %f\n", theta, phi);

    printf("init: %f     %f     %f     %f     %f     %f\n\n", rk[0], rk[1], rk[2], r_block[0], r_block[1], r_block[2]);

    select_action(epsilon, Q, s_alpha, s_mu, &a_alpha, &a_mu); // check this

    it = 0;
    reward = 0;
    tot_reward = 0;

    integration_trajectory(rk, vk, ak, r_block, v_block, a_block, \
                            r_diff, v_diff, a_diff, &theta, &phi, s_alpha, s_mu, W, &lift, &drag, &T, it);

    // ============================ KITE FLYING LOOP, STOPS WHEN FALL OCCURS ============================
    
    while (rk[2] > 0) {

      // EXIT IF MAX STEPS ARE REACHED

      if (it > max_steps){
	
        printf("MAX STEPS, %d, exiting\n", max_steps);
        printf("return=%f, space percurred=%f\n\n", tot_reward, sqrt(r_block[0]*r_block[0] + r_block[1]*r_block[1]) );
        fprintf(rew, "%d    %f\n", episode, tot_reward);

        // fill files
        fill_Q_file(Q_mat, Q);
        fill_Q_counter(Q_mat_counter, Q_counter);
        
        break;

      }

      integration_trajectory(rk, vk, ak, r_block, v_block, a_block,\
			     r_diff, v_diff, a_diff, &theta, &phi, s_alpha, s_mu, W, &lift, &drag, &T, it);
      
      reward = fabs(v_block[0]*v_block[0] + v_block[1]*v_block[1])*reward_dt;

      if (it%decision_time == 0){
	      tot_reward += reward;
      }

      // SOME CHECKS

      if (check(s_alpha, a_alpha, s_mu, a_mu) == 1) {
        printf("CHECK: s_alpha %d, s_mu %d\n", s_alpha, s_mu);
        printf("Matrix index error!!\n");
        episode = learning_episodes;
        break;
      }

      if (a_alpha > 2 || a_alpha < 0 || a_mu > 2 || a_mu < 0){
        printf("a_alpha:%d, a_mu:%d\n", a_alpha, a_mu);
        printf("ERROR IN UPDATE STATE");
        episode = learning_episodes;
        break;
      }

      if ( (episode == learning_episodes - 1) && (it%decision_time == 0) ){
	      fprintf(out, "%d       %f       %f      %f      %f      %f      %f     %f      %f      %f      %f     %f\n", it, \
		      rk[0], rk[1], rk[2], r_block[0], r_block[1], r_block[2], W[0], W[1], W[2], v_block[0], v_block[1]);
      }

      // BREAK EPISODE IF KITE FALLS OR THE BLOCK TAKES FLIGHT, UPDATING Q WITH PENALTY

      if (rk[2] <= 0. || !(m_block*g >= T*cos(theta)) ) {

        if (rk[2] <= 0.){
          printf("Kite Fall, steps %d, z<0, break\n", it);
        } else { 
          printf("Block takes flight, steps %d, break\n", it);
        }

        printf("return=%f, space percurred=%f\n\n", tot_reward, sqrt(r_block[0]*r_block[0] + r_block[1]*r_block[1]));

        Q[s_alpha*n_act*n_bank*n_act + a_alpha*n_act*n_bank + s_mu*n_act + a_mu] +=
          Alpha*(reward + penalty - Q[s_alpha*n_act*n_bank*n_act + a_alpha*n_act*n_bank + s_mu*n_act + a_mu]);
        
        Q_counter[s_alpha*n_act*n_bank*n_act + a_alpha*n_act*n_bank + s_mu*n_act + a_mu] += 1;

        // fill files
        fprintf(rew, "%d    %f\n", episode, tot_reward);
        fill_Q_file(Q_mat, Q);
        fill_Q_counter(Q_mat_counter, Q_counter);

        break;
      }

      // FIX RADIUS AND TAKE A DECISION EVERY 0.1 SEC

      if (it%decision_time == 0){

        r_diff_modulo = sqrt(r_diff[0]*r_diff[0] + r_diff[1]*r_diff[1] + r_diff[2]*r_diff[2]);

        rk1[0] = r_block[0] + (rk[0] - r_block[0])/fabs(r_diff_modulo)*R;
        rk1[1] = r_block[1] + (rk[1] - r_block[1])/fabs(r_diff_modulo)*R;
        rk1[2] = r_block[2] + (rk[2] - r_block[2])/fabs(r_diff_modulo)*R;

        rk[0] = rk1[0];
        rk[1] = rk1[1];
        rk[2] = rk1[2];
            
        // FIND STATE S1 USING ACTION A

        s_alpha1 = update_state(s_alpha, a_alpha);
        s_mu1 = update_state(s_mu, a_mu);

        // SEARCH FOR NEXT ACTION A1

        select_action(epsilon, Q, s_alpha1, s_mu1, &a_alpha1, &a_mu1);

        // UPDATE Q MATRIX

        if (it == max_steps){
          Q[s_alpha*n_act*n_bank*n_act + a_alpha*n_act*n_bank + s_mu*n_act + a_mu] += \
            Alpha*(reward - Q[s_alpha*n_act*n_bank*n_act + a_alpha*n_act*n_bank + s_mu*n_act + a_mu]);
                          
          Q_counter[s_alpha*n_act*n_bank*n_act + a_alpha*n_act*n_bank + s_mu*n_act + a_mu] += 1;
        }
        else {
          Q[s_alpha*n_act*n_bank*n_act + a_alpha*n_act*n_bank + s_mu*n_act + a_mu] += \
            Alpha*(reward + Gamma*Q[s_alpha1*n_act*n_bank*n_act + a_alpha1*n_act*n_bank + s_mu1*n_act + a_mu1] \
            - Q[s_alpha*n_act*n_bank*n_act + a_alpha*n_act*n_bank + s_mu*n_act + a_mu]);
                          
          Q_counter[s_alpha*n_act*n_bank*n_act + a_alpha*n_act*n_bank + s_mu*n_act + a_mu] += 1;
        }

        // MOVE ON: S = S1, A = A1

        s_alpha = s_alpha1;
        s_mu = s_mu1;

        a_alpha = a_alpha1;
        a_mu = a_mu1;

      }

      it += 1;
    }

    episode += 1;

  }

  save_matrix(Q, "Q_matrix.dat");

  free(rk);
  free(rk1);
  free(vk);
  free(ak);

  free(r_diff);
  free(v_diff);
  free(a_diff);

  free(Q);
  free(Q_counter);

  free(r_block);
  free(v_block);
  free(a_block);

  fclose(out);
  fclose(rew);
  fclose(Q_mat);
  fclose(Q_mat_counter);
  fclose(policy);

  return 0;

    }

#include <stdio.h>
#include <math.h>

#ifndef __sarsa__
#define __sarsa__

#define PI 3.1415926535897932384626433

#define num_saved_matrices 500
#define learning_episodes 20000

#define dim 2
#define PENALTY -200

#define n_actions 3 // 0 diminuisco alpha, 1 rimango, 2 aumento
#define n_thetas 7
#define n_thetas_vel 13

#define decision_time 1000

#define Gamma 0.9999999999
<<<<<<< HEAD
<<<<<<< HEAD
#define max_steps 2000000
=======
#define learning_episodes 3000
#define max_steps 1000000
>>>>>>> 8a4af4f473d680572248ff2f1a7b5406024ab1e4
=======
#define learning_episodes 3000
#define max_steps 1000000
>>>>>>> 8a4af4f473d680572248ff2f1a7b5406024ab1e4

#define s_alpha0 10

double Alpha = 0.001; // 10e-5 come punto di partenza
double epsilon = 0.9;

//#define DEBUGSARSA

double thetas[n_thetas] = {0., PI/8., PI/4., PI*3/8., PI/2., PI*3/4, PI};
double thetas_vel[n_thetas_vel] = {-20., -10., -2., -1.5, -1., -0.5, 0., 0.5, 1., 1.5,  2., 10., 20.};

void fill_Q_mat(FILE *Q_mat_file, double *Q, int episode){

  for (int kk=0; kk<n_thetas; kk++){
    for (int p=0; p<n_thetas_vel; p++){
      for (int i=0; i<n_alphas; i++){
        fprintf(Q_mat_file,"%d,", episode); //print episode
        fprintf(Q_mat_file,"%d,", i); // print alpha idx
        fprintf(Q_mat_file,"%f,", thetas[kk]); // print theta
        fprintf(Q_mat_file,"%f,", thetas_vel[p]); // print vtheta idx
        for (int j=0; j<n_actions; j++){
          fprintf(Q_mat_file,"%f", Q[i*n_actions*n_thetas*n_thetas_vel + j*n_thetas*n_thetas_vel + 
                                    kk*n_thetas_vel + p]);
          if (j!=n_actions-1){ fprintf(Q_mat_file,","); }
        }
        fprintf(Q_mat_file,"\n");
      }
      fprintf(Q_mat_file,"\n");
    }
    fprintf(Q_mat_file,"\n");
  }
  fprintf(Q_mat_file,"\n");
  
}

void fill_Q_count(FILE *Q_mat_file, int *Q, int episode){

  for (int kk=0; kk<n_thetas; kk++){
    for (int p=0; p<n_thetas_vel; p++){
      for (int i=0; i<n_alphas; i++){
        fprintf(Q_mat_file,"%d,", episode); //print episode
        fprintf(Q_mat_file,"%d,", i); // print alpha idx
        fprintf(Q_mat_file,"%f,", thetas[kk]); // print theta
        fprintf(Q_mat_file,"%f,", thetas_vel[p]); // print vtheta idx
        for (int j=0; j<n_actions; j++){
          fprintf(Q_mat_file,"%d", Q[i*n_actions*n_thetas*n_thetas_vel + j*n_thetas*n_thetas_vel + 
                                    kk*n_thetas_vel + p]);
          if (j!=n_actions-1){ fprintf(Q_mat_file,","); }
        }
        fprintf(Q_mat_file,"\n");
      }
      fprintf(Q_mat_file,"\n");
    }
    fprintf(Q_mat_file,"\n");
  }
  fprintf(Q_mat_file,"\n");
  
}

int check(int s_alpha, int a_alpha, int s_theta, int s_vtheta){

  #ifdef DEBUGSARSA
  printf("s_alpha %d\n", s_alpha);
  printf("s_theta %d\n", s_theta);
  printf("s_vel_theta %d\n", s_vel_theta);
  #endif

  int idx_matrix = s_alpha*n_actions*n_thetas*n_thetas_vel + a_alpha*n_thetas*n_thetas_vel +
      s_theta*n_thetas_vel + s_vtheta;
  int c = 0;
  if ( idx_matrix > n_alphas*n_actions*n_thetas*n_thetas_vel || idx_matrix < 0 ){
    c = 1;
  }
  return c;
}

void initialize_Q(double *Q, double max_Q_value){

  for (int i=0; i<n_alphas; i++){
    for (int j=0; j<n_actions; j++){
      for (int p=0; p<n_thetas; p++){
        for (int kk=0; kk<n_thetas_vel; kk++){

          Q[i*n_actions*n_thetas*n_thetas_vel + j*n_thetas*n_thetas_vel + p*n_thetas_vel + kk] = max_Q_value;

          if (i == 0 && j ==0){
              Q[i*n_actions*n_thetas*n_thetas_vel + j*n_thetas*n_thetas_vel + p*n_thetas_vel + kk] = 0;
          }
          if (i == n_alphas-1 && j == n_actions-1){
              Q[i*n_actions*n_thetas*n_thetas_vel + j*n_thetas*n_thetas_vel + p*n_thetas_vel + kk] = 0;
          }

        }
      }
    }
  }
}

void initialize_Q_count(int *Q, int var){

  for (int i=0; i<n_alphas; i++){
    for (int j=0; j<n_actions; j++){
      for (int p=0; p<n_thetas; p++){
        for (int kk=0; kk<n_thetas_vel; kk++){
          Q[i*n_actions*n_thetas*n_thetas_vel + j*n_thetas*n_thetas_vel + p*n_thetas_vel + kk] = var;
        }
      }
    }
  }
}


int find_state(double v, double * v_array){

  int slice = 0;
  for (int i=0; i<sizeof(v_array)/sizeof(v_array[0]); i++){
      if (v >= v_array[i] && v < v_array[i+1]){
          slice = i;
      }
  }

  #ifdef DEBUGSARSA
    printf("slice = %d\n", slice);
  #endif

  return slice;
}

int update_state(int s, int a){
  int s1;
  if (a == 0){
    s1 = s - 1;
  }
  else if (a == 1){
    s1 = s;
  }
  else if (a == 2){
    s1 = s + 1;
  }
  else {
    printf(" ===========================> UPDATE STATE ERROR!!!!!!!!");
    s1 = 100;
  }
  
  return s1;
}

int select_alpha_action(double epsilon, double *Q, int s_alpha, int s_theta, int s_vtheta, int episode){

  int a1;

  #ifdef DEBUGSARSA
  printf("***select mu action from s_mu = %d, s_z = %d, s_alpha = %d\n", s_alpha, s_theta, s_vel_theta);
  #endif

  double probab = (double)rand()/RAND_MAX; //rand() % 100;

  if (probab >= epsilon){ // scelgo la nuova azione come numero random

    #ifdef DEBUGSARSA
        printf("   random\n");
    #endif
    if (s_alpha == 0){
        a1 = (rand() % 2) + 1;
    }
    else if (s_alpha == n_alphas-1){
        a1 = (rand() % 2); // random number da 0 a 2
    }
    else {
        a1 = rand() % 3;
    }
  }

  else {

    a1 = 0;
    if (s_alpha == 0) {
      a1 = 1;
    }

    int max_value = n_actions;
    if (s_alpha == n_alphas - 1) { 
      max_value = n_actions - 1; 
    }

    #ifdef DEBUGSARSA
        printf("  starting action a1_0 = %d, max value = %d\n", a1, max_value);
    #endif

    /* LOOP TO SEARCH NEXT ACTION */ 

    for (int i=a1+1; i<max_value; i++){
      if ( Q[s_alpha*n_actions*n_thetas*n_thetas_vel + i*n_thetas*n_thetas_vel + s_theta*n_thetas_vel + s_vtheta]> 
             Q[s_alpha*n_actions*n_thetas*n_thetas_vel + a1*n_thetas*n_thetas_vel + s_theta*n_thetas_vel + s_vtheta]) { 
        a1 = i;
      }
    }
  }

  #ifdef DEBUGSARSA
  printf("  a1 = %d\n**end select action\n", a1);
  #endif

  return a1;
}



#endif

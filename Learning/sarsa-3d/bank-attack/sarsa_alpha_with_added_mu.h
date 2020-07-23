#include <stdio.h>
#include <math.h>

#ifndef __sarsa__
#define __sarsa__

#define PI 3.1415926535897932384626433

#define num_saved_matrices 500
#define learning_episodes 500000

#define PENALTY -300.0

#define n_actions 3 // 0 diminuisco alpha, 1 rimango, 2 aumento
#define n_bank 13

#define decision_time 1000

#define Gamma 0.9999999999
#define max_steps 10000000

#define s_alpha0 12

#ifndef n_alphas
#define n_alphas 15
#endif

// n bank 13, -10, -8, -4, -6, -2, -1, 0, 1, 2, 4, 6, 8, 10
double bank[n_bank] = {-0.174532925, -0.13962634016, -0.10471975512, -0.06981317008, -0.03490658504, -0.01745329252, 0., 0.01745329252, 0.03490658504, 0.06981317008, 0.10471975512, 0.13962634016, 0.174532925};
#define s_bank0 6

double Alpha = 0.5; // ordine di grandezza=h/tempo_episodio_in_sec(al minimo caduta)
double epsilon = 0.05;

void fill_Q_mat(FILE *Q_mat_file, double *Q, int episode){

  for (int i=0; i<n_alphas; i++){
  for (int j=0; j<n_bank; j++){
    for (int ii=0; ii<n_actions; ii++){
    for (int jj=0; jj<n_actions; jj++){
    fprintf(Q_mat_file,"%d,", episode);
    fprintf(Q_mat_file,"%d,", i);
    fprintf(Q_mat_file,"%d,", j);
    fprintf(Q_mat_file,"%d,", ii);
    fprintf(Q_mat_file,"%d,", jj);
      fprintf(Q_mat_file,"%f\n", Q[j*n_actions*n_actions*n_alphas + jj*n_actions*n_alphas + i*n_actions + ii]);
    }
    }
  }
  }
  
}

void fill_Q_count(FILE *Q_mat_file, int *Q, int episode){

  for (int i=0; i<n_alphas; i++){
  for (int j=0; j<n_bank; j++){
    for (int ii=0; ii<n_actions; ii++){
    for (int jj=0; jj<n_actions; jj++){
    fprintf(Q_mat_file,"%d,", episode);
    fprintf(Q_mat_file,"%d,", i);
    fprintf(Q_mat_file,"%d,", j);
    fprintf(Q_mat_file,"%d,", ii);
    fprintf(Q_mat_file,"%d,", jj);
      fprintf(Q_mat_file,"%d\n", Q[j*n_actions*n_actions*n_alphas + jj*n_actions*n_alphas + i*n_actions + ii]);
    }
    }
  }
  }
  
}

void load_matrix_3d(double * Q, char * name){

  FILE *dat = fopen(name, "r"); // opens file for reading
  if (dat)
  {
    for (int i=0; i<n_alphas; i++){
    for (int j=0; j<n_bank; j++){
      for (int ii=0; ii<n_actions; ii++){
      for (int jj=0; jj<n_actions; jj++){
        fscanf(dat, "%lf ", &Q[j*n_actions*n_actions*n_alphas + jj*n_alphas*n_actions + i*n_actions + ii]);
      }
      }
    }
    }

  } 
  else { printf("error loading file\n"); }

  fclose(dat);

}

void load_matrix_only_alpha(double * Q, char * name){

  int idx_mu = 1;

  FILE *dat = fopen(name, "r"); // opens file for reading
  if (dat)
  {
    for (int i = 0; i < n_alphas; i++){
      for (int j = 0; j < n_actions; j++){
        // non ho messo 0 perche` per mu=0 e azione=0 la matrice e` tutta nulla!
        fscanf(dat, "%lf ", &Q[idx_mu*n_actions*n_actions*n_alphas + 0*n_alphas*n_actions + i*n_actions + j]);
      }
    }

    int kk = 0;
    int ll = 0;
    while(kk < n_bank){
      ll = 0;
    while (ll < n_actions){

      for (int i=0; i<n_alphas; i++){
        for (int j=0; j<n_actions; j++){
          if (kk == 0 && ll==0){
            Q[kk*n_alphas*n_actions*n_actions + ll*n_alphas*n_actions + i*n_actions + j]  = 0;
          }
          else if (kk == n_bank-1 && ll == n_actions-1){
            Q[kk*n_alphas*n_actions*n_actions + ll*n_alphas*n_actions + i*n_actions + j]  = 0;
          } else {
            Q[kk*n_alphas*n_actions*n_actions + ll*n_alphas*n_actions + i*n_actions + j] = \
            Q[idx_mu*n_actions*n_actions*n_alphas + 0*n_alphas*n_actions + i*n_actions + j];
          }
          }
        }
        ll++;
      }
    kk++;
    }
  }
  else { printf("error loading file\n"); }

  fclose(dat);
}

void print_mat(double * Q){

  for (int j=0; j<n_bank; j++){
    for (int jj=0; jj<n_actions; jj++){
    printf("bank=%f, action=%d\n", bank[j], jj);
    for (int i=0; i<n_alphas; i++){
      for (int ii=0; ii<n_actions; ii++){
      printf("%f ", Q[j*n_actions*n_actions*n_alphas + jj*n_actions*n_alphas + i*n_actions + ii]);
    }
    printf("\n");
    }
    printf("\n\n");
  }
  }
  printf("\n\n");
  
}
/*
void save_matrix(double * Q, char * name){

  FILE *dat = fopen(name, "w"); // opens new file for writing
  if (dat)
  {
    for (int kk = 0; kk < n_alphas; kk++){
      for (int j = 0; j < n_actions; j++){
        fprintf(dat, "%16.8e ", Q[kk*n_actions + j]);
      }
    }
  }

  fclose(dat);
}

void load_matrix(double * Q, char * name){

  FILE *dat = fopen(name, "r"); // opens file for reading //"Qmatrix.dat"
  if (dat)
  {
    for (int kk = 0; kk < n_alphas; kk++){
      for (int j=0; j<n_actions; j++){
        fscanf(dat, "%lf", &Q[kk*n_actions + j] );
      }
    }
  }
  else { printf("error loading file\n"); }

  fclose(dat);
}*/

int check(int s_alpha, int s_bank, int a_alpha, int a_bank){

    #ifdef DEBUGSARSA
    //printf("s_alpha %d\n", s_alpha);
    #endif
    int c = 0;
    int idx = s_bank*n_actions*n_actions*n_alphas + a_bank*n_actions*n_alphas + s_alpha*n_actions + a_alpha;
    if ( (idx >= n_actions*n_actions*n_alphas*n_bank || idx < 0)){
      c = 1;
    }
    return c;
}

void initialize_Q(double * Q, double max_Q_value){

    for (int i=0; i<n_alphas; i++){
    for (int j=0; j<n_bank; j++){
      for (int ii=0; ii<n_actions; ii++){
      for (int jj=0; jj<n_actions; jj++){
        if (i == 0 && ii == 0) Q[j*n_actions*n_actions*n_alphas + jj*n_alphas*n_actions + i*n_actions + ii] = 0;
        if (j == 0 && jj == 0) Q[j*n_actions*n_actions*n_alphas + jj*n_alphas*n_actions + i*n_actions + ii] = 0;
        if (i == n_alphas-1 && ii == n_actions-1) Q[j*n_actions*n_actions*n_alphas + jj*n_alphas*n_actions + i*n_actions + ii] = 0;
        if (j == n_bank-1 && jj == n_actions-1) Q[j*n_actions*n_actions*n_alphas + jj*n_alphas*n_actions + i*n_actions + ii] = 0;
        else Q[j*n_actions*n_actions*n_alphas + jj*n_alphas*n_actions + i*n_actions + ii] = max_Q_value;
      }
      }
    }
    }
}

void initialize_Q_count(int * Q, int value){

    for (int i=0; i<n_bank*n_actions*n_alphas*n_actions; i++){
      Q[i] = 0;
    }
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

void select_action(double epsilon, double * Q, int s_alpha, int *a_alpha, int s_bank, int *a_bank){

  double probab = (double)rand()/RAND_MAX;
  
  if (probab <= epsilon){ // random action
  
    if (s_alpha == 0){
        *a_alpha = (rand() % 2) + 1; // random number between 1 and 2 (compresi)
    }
    else if (s_alpha == n_alphas-1){
        *a_alpha = (rand() % 2); // random number da 0 a 2
    }
    else {
        *a_alpha = rand() % 3;
    }

    if (s_bank == 0){
      *a_bank = (rand() % 2) + 1; // random number between 1 and 2 (compresi)
    }
    else if (s_bank == n_bank-1){
      *a_bank = (rand() % 2); // random number da 0 a 2
    }
    else {
      *a_bank = rand() % 3;
    }

    }  else { // GREEDY ACTION
    //printf("not random\n");

    /* fix starting action for the search: 0 in any case, and 1 only if the 
      starting state is the beginning of the array */

    int starting_action_alpha = 1;
    int starting_action_mu = 1;
    *a_bank = starting_action_alpha;
    *a_alpha = starting_action_mu;
  
    int a_alpha_min = 0;
    int a_bank_min = 0;

    if (s_alpha == 0) {
      a_alpha_min = 1;
    }
    if (s_bank == 0) {
      a_bank_min = 1; 
    }

    /* fix ending state: 3 (n_actions), NOT INCLUDED, in any state; 2 (n_actions -1)
     only if starting state is the last one in the array */

    int a_alpha_max = n_actions;
    int a_bank_max = n_actions;

    if (s_alpha == n_alphas - 1) { 
      a_alpha_max = n_actions - 1; 
    }
    if (s_bank == n_bank - 1) { 
      a_bank_max = n_actions - 1; 
    }

    /* LOOP TO SEARCH NEXT ACTION */ 
    //printf("a_alpha_min = %d, a_alpha_max=%d, a_bank_min=%d, a_bank_max=%d\n", a_alpha_min, a_alpha_max, a_bank_min, a_bank_max);
    
    for (int ii=a_alpha_min; ii<a_alpha_max; ii++){
      for (int jj=a_bank_min; jj<a_bank_max; jj++){
     
        if (Q[s_bank*n_actions*n_alphas*n_actions + jj*n_actions*n_alphas + s_alpha*n_actions + ii] >
            Q[s_bank*n_actions*n_alphas*n_actions + *a_bank*n_actions*n_alphas + s_alpha*n_actions + *a_alpha])  { 
          *a_alpha = ii;
          *a_bank = jj;
          //printf("Modify! New alpha action=%d, new mu actions=%d\n", *a_alpha1, *a_mu1);
        }
 
      }
      }
    }

    //printf("chosen: a_alpha=%d, a_bank=%d\n", *a_alpha, *a_bank);


}

double expected_update(double * Q, int s_alpha, int a_alpha, int s_bank, int a_bank){

  double update = 0;

  return update;
}

#endif
#include <stdio.h>
#include <math.h>

#ifndef __sarsa__
#define __sarsa__

#define PI 3.1415926535897932384626433

#define num_saved_matrices 400
#define learning_episodes 150000

#define PENALTY -300.0

#define n_actions 3 // 0 diminuisco alpha, 1 rimango, 2 aumento
#define n_angles 30 //18

#define decision_time 1000

#define Gamma 0.9999999999
#define max_steps 2000000

#define s_alpha0 10

#ifndef n_alphas
#define n_alphas 15
#endif

double Alpha = 0.5; // ordine di grandezza=h/tempo_episodio_in_sec(al minimo caduta)
double epsilon = 0.05;

double angles[n_angles+1] = {-4, -3.5, -3.2, -3.1, -3.0, -2.9, -2.7, -2.4, -2.1, -1.8, -1.5, -1.2, -0.9, -0.6, \
-0.3, 0., 0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.1, 2.4, 2.7, 2.9, 3.0, 3.1, 3.2, 3.5, 4.}; // impose value n_angles to 29/// I imposed it to 30!!

//#define DEBUGSARSA

void fill_Q_mat(FILE *Q_mat_file, double *Q, int episode){

  for (int kk=0; kk<n_angles; kk++){
  for (int i=0; i<n_alphas; i++){
    fprintf(Q_mat_file,"%d,", episode);
    fprintf(Q_mat_file,"%d,", kk);
    fprintf(Q_mat_file,"%d,", i);
    for (int j=0; j<n_actions; j++){
      fprintf(Q_mat_file,"%f", Q[kk*n_alphas*n_actions + i*n_actions + j]);
      if (j!=n_actions-1){ 
        fprintf(Q_mat_file,","); 
      }
    }
  fprintf(Q_mat_file,"\n");
  }
  fprintf(Q_mat_file,"\n");
  }
  
}

void fill_Q_count(FILE *Q_mat_file, int *Q, int episode){

  for (int kk=0; kk<n_angles; kk++){
  for (int i=0; i<n_alphas; i++){
    fprintf(Q_mat_file,"%d,", episode);
    fprintf(Q_mat_file,"%d,", kk);
    fprintf(Q_mat_file,"%d,", i);
    for (int j=0; j<n_actions; j++){
      fprintf(Q_mat_file,"%d", Q[kk*n_alphas*n_actions + i*n_actions + j]);
      if (j!=n_actions-1){ fprintf(Q_mat_file,","); }
    }
  fprintf(Q_mat_file,"\n");
  }
  fprintf(Q_mat_file,"\n");
  }
  
}

void save_matrix(double * Q, char * name){

  FILE *dat = fopen(name, "w"); // opens new file for writing
  if (dat)
  {
    for (int kk = 0; kk < n_angles; kk++){
    for (int i = 0; i < n_alphas; i++){
      for (int j = 0; j < n_actions; j++){
        fprintf(dat, "%16.8e ", Q[kk*n_alphas*n_actions + i*n_actions + j]);
      }
    }
    }
  }

  fclose(dat);
}

void load_matrix(double * Q, char * name){

  FILE *dat = fopen(name, "r"); // opens file for reading
  if (dat)
  {
    for (int kk = 0; kk < n_angles; kk++){
      for (int i = 0; i < n_alphas; i++){
        for (int j = 0; j < n_actions; j++){
          fscanf(dat, "%lf ", &Q[kk*n_alphas*n_actions + i*n_actions + j]);
        }
      }
    }
  }
  else { printf("Nic: error loading file\n"); }

  fclose(dat);
}

void load_counter(int * Q, char * name){

  FILE *dat = fopen(name, "r"); // opens file for reading
  if (dat)
  {
    for (int kk = 0; kk < n_angles; kk++){
      for (int i = 0; i < n_alphas; i++){
        for (int j = 0; j < n_actions; j++){
          fscanf(dat, "%d ", &Q[kk*n_alphas*n_actions + i*n_actions + j]);
        }
      }
    }
  }
  else { printf("Nic: error loading file\n"); }

  fclose(dat);
}

int check(int s_alpha, int action, int s_vrel_mod){

    #ifdef DEBUGSARSA
    //printf("s_alpha %d\n", s_alpha);
    #endif

    int idx = s_vrel_mod*n_alphas*n_actions + s_alpha*n_actions + action;
    int c = 0;
    if ( (idx >= n_alphas*n_actions*n_angles || idx < 0)){
      c = 1;
    }
    return c;
}

void initialize_Q(double * Q, double max_Q_value){

  for (int kk = 0; kk<n_angles; kk++){
    for (int i=0; i<n_alphas; i++){
      for (int j=0; j<n_actions; j++){
        if (i == 0 && j == 0) Q[kk*n_alphas*n_actions + i*n_actions + j] = 0;
        else if (i == n_alphas-1 && j == n_actions-1) Q[kk*n_alphas*n_actions + i*n_actions + j] = 0;
        else Q[kk*n_alphas*n_actions + i*n_actions + j] = max_Q_value;
      }
    }
  }
}


void initialize_pi(double * pi){

  for (int kk = 0; kk<n_angles; kk++){
    for (int i=0; i<n_alphas; i++){
      for (int j=0; j<n_actions; j++){
        if (i == 0){
          if (j == 0) { pi[kk*n_alphas*n_actions + i*n_actions + j] = 0; }
          else { pi[kk*n_alphas*n_actions + i*n_actions + j] = 0.5; }
        }
        else if (i == n_alphas - 1){
          if (j == n_actions - 1){ pi[kk*n_alphas*n_actions + i*n_actions + j] = 0; }
          else { pi[kk*n_alphas*n_actions + i*n_actions + j] = 0.5; }
        }
        else pi[kk*n_alphas*n_actions + i*n_actions + j] = 0.3333333333;
      }
    }
  }
}

void initialize_Q_count(int * Q, int value){

    for (int i=0; i<n_angles*n_actions*n_alphas; i++){
      Q[i] = 0;
    }
}

int find_state_angle(double v){

  int slice = 0;

  for (int i=0; i<n_angles; i++){
      if (v >= angles[i] && v < angles[i+1]){
          slice = i;
      }
  }

  #ifdef DEBUGSARSA
    printf("slice = %d\n", slice);
  #endif

  return slice;
}

void print_matrix(double * mat){

  for (int kk=0; kk<n_angles; kk++){
    printf("angle between %f and %f\n", angles[kk], angles[kk+1]);
    for (int i=0; i<n_alphas; i++){
      printf("alph idx=%d ", i);
      for (int j=0; j<n_actions; j++){
        printf("%f ",mat[kk*n_alphas*n_actions + i*n_actions + j] );
      }
      printf("\n");
    }
    printf("\n");
  }

}

void print_counter(int * Q){

  for (int kk=0; kk<n_angles; kk++){
    printf("angle between %f and %f\n", angles[kk], angles[kk+1]);
    for (int i=0; i<n_alphas; i++){
      printf("alph idx=%d ", i);
      for (int j=0; j<n_actions; j++){
        printf("%d ", Q[kk*n_alphas*n_actions + i*n_actions + j] );
      }
      printf("\n");
    }
    printf("\n");
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

int select_alpha_action(double epsilon, double * Q, int s_alpha,int s_vrel_mod){

  int a1;

  double probab = (double)rand()/RAND_MAX;
  
  if (probab <= epsilon){ // random action
  
    if (s_alpha == 0){
        a1 = (rand() % 2) + 1; // random number between 1 and 2 (compresi)
    }
    else if (s_alpha == n_alphas-1){
        a1 = (rand() % 2); // random number da 0 a 2
    }
    else {
        a1 = rand() % 3;
    }

  }  else { // GREEDY ACTION

    /* fix starting action for the search: 0 in any case, and 1 only if the 
      starting state is the beginning of the array */

    int starting_value = 1;
    a1 = starting_value;

    int min_value = 0;
    if (s_alpha == 0) { 
      min_value = 1;
    }

    /* fix ending state: 3 (n_actions), NOT INCLUDED, in any state; 2 (n_actions -1)
     only if starting state is the last one in the array */

    int max_value = n_actions;
    if (s_alpha == n_alphas - 1) { 
      max_value = n_actions - 1; 
    }

    /* LOOP TO SEARCH NEXT ACTION */ 

    for (int i=min_value; i<max_value; i++){
      #ifdef DEBUGSARSA
      printf("   Q[(s_mu+%d-1)] = %f\n", i, Q[(s_alpha+i-1)]);
      #endif
      // io sono in Q[s_alpha*n_actions] 
      if (Q[s_vrel_mod*n_alphas*n_actions + s_alpha*n_actions + i] > Q[s_vrel_mod*n_alphas*n_actions + s_alpha*n_actions + a1]) { 
        a1 = i;
      }
    }
  }


  return a1;
}

double expected_update(double * Q, double * pi, int s_alpha1, int s_vrel_angle1){

  int min_action = 0;
  int max_action = 3;
  if (s_alpha1 == 0){
    min_action = 1;
  }
  if (s_alpha1 == n_alphas - 1){
    max_action = n_actions - 1;
  }
    
  double update = 0;
  for (int i=min_action; i<max_action; i++){
    update += Q[s_vrel_angle1*n_alphas*n_actions + s_alpha1*n_actions + i]*pi[s_vrel_angle1*n_alphas*n_actions + s_alpha1*n_actions + i];
  }

  return update;

}


void softmax_policy(double *pi, double *Q, int *Q_count){

double sum = 0;
for (int kk=0; kk<n_angles; kk++){ // numero di angoli della matrice
  for (int i=0; i<n_alphas; i++){ // numero di attack 
    
    for (int j=0; j<n_actions; j++){ // numero di azioni

    sum = 0;
    if ( i == 0 ){
      for (int jj=1; jj<n_actions; jj++){
        sum += exp(Q[kk*n_alphas*n_actions + i*n_actions + j] - Q[kk*n_alphas*n_actions + i*n_actions + jj]);
      }
      pi[kk*n_alphas*n_actions + i*n_actions + j] = 1/sum;
      pi[kk*n_alphas*n_actions + i*n_actions + 0] = 0;

    }

    else if ( i == n_alphas-1 ){
      pi[kk*n_alphas*n_actions + i*n_actions + j] = 0;
      for (int jj=0; jj<n_actions-1; jj++){
        sum += exp(Q[kk*n_alphas*n_actions + i*n_actions + j] - Q[kk*n_alphas*n_actions + i*n_actions + jj]);
      }
      pi[kk*n_alphas*n_actions + i*n_actions + j] = 1/sum;
      pi[kk*n_alphas*n_actions + i*n_actions + n_actions-1] = 0;
    }

    else{
      for (int jj=0; jj<n_actions; jj++){
        sum += exp(Q[kk*n_alphas*n_actions + i*n_actions + j] - Q[kk*n_alphas*n_actions + i*n_actions + jj]);
      }
      pi[kk*n_alphas*n_actions + i*n_actions + j] = 1/sum;
    }
    }
  }
}

}


int softmax_check(double *pi){

  int checker = 0;

  double sum = 0;
  for (int kk=0; kk<n_angles; kk++){ // numero di angoli della matrice
    for (int i=0; i<n_alphas; i++){

      sum = 0;
      
      for (int j=0; j<n_actions; j++){ 
          sum += pi[kk*n_alphas*n_actions + i*n_actions + j];
      }
      if ( fabs(sum - 1) > 0.1){
        printf("sum=%f\n", sum);
        printf("fabs(sum - 1)=%f\n", fabs(sum - 1));
        checker = 1;
      }
    }
  }

  return checker;
}


#endif
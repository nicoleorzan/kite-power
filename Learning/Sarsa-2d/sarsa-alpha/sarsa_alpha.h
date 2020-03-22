#include <stdio.h>
#include <math.h>

#ifndef __sarsa__
#define __sarsa__

#define PI 3.1415926535897932384626433

#define dim 2
#define PENALTY -200

#define n_actions 3 // 0 diminuisco alpha, 1 rimango, 2 aumento

#define decision_time 1000

#define Gamma 0.9999999999
#define learning_episodes 3000
#define max_steps 1000000

#define s_alpha0 10

double Alpha = 0.001; // ordine di grandezza=h/tempo_episodio_in_sec(al minimo caduta)
double epsilon = 0.9;

//n_alphas defined in fagiano_model_constant

//#define DEBUGSARSA

void fill_Q_mat(FILE *Q_mat_file, double *Q, int episode){

  for (int i=0; i<n_alphas; i++){
    fprintf(Q_mat_file,"%d,", episode);
    fprintf(Q_mat_file,"%d,", i);
    for (int j=0; j<n_actions; j++){
      fprintf(Q_mat_file,"%f", Q[i*n_actions + j]);
      if (j!=n_actions-1){ fprintf(Q_mat_file,","); }
    }
  fprintf(Q_mat_file,"\n");
  }
  fprintf(Q_mat_file,"\n");
  
}

void fill_Q_count(FILE *Q_mat_file, int *Q){

  for (int i=0; i<n_alphas; i++){
    for (int j=0; j<n_actions; j++){
          fprintf(Q_mat_file,"%d ", Q[i*n_actions + j]);
        }
      }
  fprintf(Q_mat_file,"\n");
  
}

void print_mat(double * Q){

  for (int i=0; i<n_alphas; i++){
    printf("alpha=%.1f  ", alphas[i]);
    for (int j=0; j<n_actions; j++){
      printf("%f ", Q[i*n_actions + j]);
    }
    printf("\n");
  }
  printf("\n\n");
  
}

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
}

int check(int s_alpha){

    #ifdef DEBUGSARSA
    //printf("s_alpha %d\n", s_alpha);
    #endif

    int c = 0;
    if ( (s_alpha >= n_alphas*n_actions || s_alpha < 0)){
      c = 1;
    }
    return c;
}

void initialize_Q(double * Q, double max_Q_value){

    for (int i=0; i<n_alphas; i++){
      for (int j=0; j<n_actions; j++){
        if (i == 0 && j == 0) Q[i*n_actions + j] = 0;
        else if (i == n_alphas-1 && j == n_actions-1) Q[i*n_actions + j] = 0;
        else Q[i*n_actions + j] = max_Q_value;
      }
    }
}

void initialize_Q_count(int * Q, int value){

    for (int i=0; i<n_alphas; i++){
      for (int j=0; j<n_actions; j++){
        Q[i*n_actions + j] = value;
      }
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

int select_alpha_action(double epsilon, double * Q, int s_alpha, int episode){

  int a1;

  #ifdef DEBUGSARSA
  printf("\n***select action from status = %d\n", s_alpha);
  #endif

  double probab = (double)rand()/RAND_MAX;
  
  if (probab >= epsilon){ // new action as a random number
   
    #ifdef DEBUGSARSA
    if (episode == learning_episodes-1){
        printf("    random\n");
        printf("    probab= %f, epsilon=%f\n", probab, epsilon);
    }
    #endif
    if (s_alpha == 0){
        a1 = (rand() % 2) + 1; // random number between 1 and 2 (compresi)
    }
    else if (s_alpha == n_alphas-1){
        a1 = (rand() % 2); // random number da 0 a 2
    }
    else {
        a1 = rand() % 3;
    }
    }  else {
    
    #ifdef DEBUGSARSA
    if (episode == learning_episodes-1){
        printf("    not random\n");
        printf("    probab= %f, epsilon=%f\n", probab, epsilon);
    }
    #endif

    /* fix starting action for the search: 0 in any case, and 1 only if the 
      starting state is the beginning of the array */

    a1 = 0;
    if (s_alpha == 0) { 
      a1 = 1; 
    }

    /* fix ending state: 3 (n_actions), NOT INCLUDED, in any state; 2 (n_actions -1)
     only if starting state is the last one in the array */

    int max_value = n_actions;
    if (s_alpha == n_alphas - 1) { 
      max_value = n_actions - 1; 
    }

    #ifdef DEBUGSARSA
        printf("  starting action a1_0 = %d, max value (NOT INCLUDED)= %d\n", a1, max_value);
    #endif

    #ifdef DEBUGSARSA
    printf("  Q[(s_alpha+%d-1)] = %f\n", a1, Q[(s_alpha+a1-1)]);
    #endif

    /* LOOP TO SEARCH NEXT ACTION */ 

    for (int i=a1+1; i<max_value; i++){
      #ifdef DEBUGSARSA
      printf("   Q[(s_mu+%d-1)] = %f\n", i, Q[(s_alpha+i-1)]);
      #endif
      // io sono in Q[s_alpha*n_actions] 
      if (Q[s_alpha*n_actions + i] > Q[s_alpha*n_actions + a1]) { 
        a1 = i;
      }
    }
  }

  #ifdef DEBUGSARSA
  printf("  a1 = %d\n**end select action\n\n", a1);
  #endif

  return a1;
}



#endif
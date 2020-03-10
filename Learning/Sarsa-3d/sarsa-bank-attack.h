#include <stdio.h>
#include <math.h>

#ifndef __sarsa__
#define __sarsa__

#define PI 3.1415926535897932384626433

#define n_act 3 // 0 diminuisco alpha, 1 rimango, 2 aumento
#define Gamma 0.9999999999
#define learning_episodes 20000
#define max_steps 30000000
#define n_bank 13

//#define DEBUGSARSA
//provare a variare alpha ed epsilon con matrice piccola per velocizzare il learning e poi trovare una
//discretizzazione ideale degli angoli theta

double Alpha = 0.1; // 0.01, 10e-5 come punto di partenza
double epsilon = 0.8;

double bank_angles[n_bank] = {-15, -12.5, -10, -7.5, -5, -2.5, 0, 2.5, 5, 7.5, 10, 12.5, 15};

/*void print_mat(double *Q){

  for (int i=0; i<n_alphas; i++){
    printf("alpha=%.1f  \n", alphas[i]);
    for (int j=0; j<n_actions; j++){
      printf("action %d  ", j);
      for (int k=0; k<n_bank; k++){
        printf("%f ", Q[i*n_actions*n_theta + j*n_theta + k]);
      }
      printf("\n");
    }
    printf("\n");
  }
  printf("\n\n");
  
}*/

void save_matrix(double *Q, char * name){

  FILE *dat = fopen(name, "w"); // opens new file for writing
  if (dat)
  {
    for (int i=0; i<n_alphas; i++){
      for (int j=0; j<n_act; j++){
        for (int k=0; k<n_bank; k++){
           for (int l=0; l<n_act; l++){
            fprintf(dat, "%16.8e ", Q[i*n_act*n_bank*n_act + j*n_bank*n_act + k*n_act + l]);
          }
        }
      }
    }
  }

  fclose(dat);
}

void load_matrix(double *Q, char * name){

  FILE *dat = fopen(name, "r"); // opens file for reading //"Qmatrix.dat"
  if (dat)
  {
    for (int i=0; i<n_alphas; i++){
      for (int j=0; j<n_act; j++){
        for (int k=0; k<n_bank; k++){
          for (int l=0; l<n_act; l++){
          fscanf(dat, "%lf", &Q[i*n_act*n_bank*n_act + j*n_bank*n_act + k*n_act + l]);
          }
        }
      }
    }
  }
  else { printf("error loading file\n"); }

  fclose(dat);
}

int check(int s_alpha, int a_alpha, int s_mu, int a_mu){

  #ifdef DEBUGSARSA
  //printf("s_alpha %d\n", s_alpha);
  #endif

  int c = 0;
  int matrix_index = s_alpha*n_act*n_bank*n_act + a_alpha*n_bank*n_act + s_mu*n_act + a_mu;
  if ( (matrix_index >= n_alphas*n_act*n_bank*n_act || matrix_index < 0)){
    printf("s_alpha*n_actions*n_theta + action*n_theta + s_theta = %d\n", matrix_index);
    c = 1;
  }
  return c;
}

void initialize_Q(double *Q, double max_Q_value){

    for (int i=0; i<n_alphas; i++){
      for (int j=0; j<n_act; j++){
        for (int k=0; k<n_bank; k++){
          for (int l=0; l<n_act; l++){
            if ( (i == 0 && j == 0) || (i==n_alphas-1 && j==n_act-1)){
              Q[i*n_act*n_bank*n_act + j*n_bank*n_act + k*n_act + l] = 0;
            }
            else Q[i*n_act*n_bank*n_act + j*n_bank*n_act + k*n_act + l] = max_Q_value;
           }
         }
       }
    }

}

void initialize_Q_counter(int *Qc, int value){

    for (int i=0; i<n_alphas*n_act*n_act*n_bank; i++){
        Qc[i] = 0;
    }
}

/*int find_state(double theta, double *arr){

  int slice = 0;
  for (int i=0; i<n_theta; i++){
      if (theta >= thetas[i] && theta < thetas[i+1]){
          slice = i;
      }
  }

  #ifdef DEBUGSARSA
  //printf("slice = %d\n", slice);
  #endif

  return slice;
}*/

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

void select_random(int s, int *a, int s_array_len){
    if (s == 0){
        *a = (rand() % 2) + 1; // random number between 1 and 2 (compresi)
    }
    else if (s == s_array_len-1){
        *a = (rand() % 2); // random number da 0 a 2
    }
    else {
        *a = rand() % 3;
    }
}

void select_action(double epsilon, double *Q, int s_alpha, int s_mu, int *a_alpha1, int *a_mu1){

  #ifdef DEBUGSARSA
  printf("\n***select action from alpha = %d, mu = %d\n", s_alpha, s_mu);
  #endif

  double probab = (double)rand()/RAND_MAX;

  if (probab >= epsilon){ // new action as a random number
    #ifdef DEBUGSARSA
        printf("   random\n");
    #endif

    if (s_alpha == 0){
        *a_alpha1 = (rand() % 2) + 1; // random number between 1 and 2 (compresi)
    }
    else if (s_alpha == n_alphas-1){
        *a_alpha1 = (rand() % 2); // random number da 0 a 2
    }
    else {
        *a_alpha1 = rand() % 3;
    }

    if (s_mu == 0){
        *a_mu1 = (rand() % 2) + 1; // random number between 1 and 2 (compresi)
    }
    else if (s_mu == n_bank-1){
        *a_mu1 = (rand() % 2); // random number da 0 a 2
    }
    else {
        *a_mu1 = rand() % 3;
    }

    //select_random(s_alpha, &a_alpha1, n_alphas);
    //select_random(s_mu, &a_mu1, n_bank);
        
  }

  else {

    /* fix starting action for the search: 0 in any case, and 1 only if the 
      starting state is the beginning of the array */

    *a_alpha1 = 0;
    *a_mu1 = 0;
    if (s_alpha == 0) { 
      *a_alpha1 = 1; 
    }
    if (s_mu == 0) { 
      *a_mu1 = 1; 
    }

    /* fix ending state: 3 (n_actions), NOT INCLUDED, in any state; 2 (n_actions -1)
     only if starting state is the last one in the array */

    int max_alpha = n_act;
    int max_mu = n_act;
    if (s_alpha == n_alphas - 1) { 
      max_alpha = n_act - 1; 
    }
    if (s_mu == n_bank - 1) { 
      max_mu = n_act - 1; 
    }

    #ifdef DEBUGSARSA
        printf("  starting action a1_alpha0 = %d, max alpha (NOT INCLUDED)= %d\n", a_alpha1, max_alpha);
	      printf("  starting action a1_mu0 = %d, max mu (NOT INCLUDED)= %d\n", a_mu1, max_mu);
    #endif

    /* LOOP TO SEARCH NEXT ACTION */ 

    for (int i=*a_alpha1+1; i<max_alpha; i++){
      for (int j=*a_mu1+1; j<max_mu; j++){
        #ifdef DEBUGSARSA
          printf("   Q[(s_alpha+%d-1)*n_theta + s_theta] = %f\n", i, Q[(s_alpha+i-1)*n_theta + s_theta]);
        #endif
        if (Q[s_alpha*n_act*n_bank*n_act + i*n_bank*n_act + s_mu*n_act + j] >
            Q[s_alpha*n_act*n_bank*n_act + *a_alpha1*n_bank*n_act + s_mu*n_act + *a_mu1])  { 
          *a_alpha1 = i;
          *a_mu1 = j;
        }
      }
    }
  }

  #ifdef DEBUGSARSA
  printf("  a1 = %d\n**end select action\n\n", a1);
  #endif

}



#endif

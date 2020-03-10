#include <stdio.h>
#include <math.h>

#ifndef __sarsa__
#define __sarsa__

#define PI 3.1415926535897932384626433

#define n_actions 3 // 0 diminuisco alpha, 1 rimango, 2 aumento
#define Gamma 0.9999999999
#define learning_episodes 5000
#define max_steps 3000000
#define n_theta 14

//#define DEBUGSARSA
//provare a variare alpha ed epsilon con matrice piccola per velocizzare il learning e poi trovare una
//discretizzazione ideale degli angoli theta

double Alpha = 0.1; // 0.01, 10e-5 come punto di partenza
double epsilon = 0.8;

double thetas[n_theta] = {0., PI/24., PI*2./24., PI*3/24., PI*4./24., PI*5./24., PI*6/24, PI*7/24, \
PI*8/24, PI*9/24, PI*10/24, PI*11/24, PI*12/24., PI};
//double thetas[n_theta] = {0., PI/12., PI*2./12., PI/4., PI*4./12., PI*5./12.,  PI/2., PI};

void print_mat(double *Q){

  for (int i=0; i<n_alphas; i++){
    printf("alpha=%.1f  \n", alphas[i]);
    for (int j=0; j<n_actions; j++){
      printf("action %d  ", j);
      for (int k=0; k<n_theta; k++){
        printf("%f ", Q[i*n_actions*n_theta + j*n_theta + k]);
      }
      printf("\n");
    }
    printf("\n");
  }
  printf("\n\n");
  
}

void save_matrix(double *Q, char * name){

  FILE *dat = fopen(name, "w"); // opens new file for writing
  if (dat)
  {
    for (int i=0; i<n_alphas; i++){
      for (int j=0; j<n_actions; j++){
        for (int k=0; k<n_theta; k++){
          fprintf(dat, "%16.8e ", Q[i*n_actions*n_theta + j*n_theta + k]);
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
      for (int j=0; j<n_actions; j++){
        for (int k=0; k<n_theta; k++){
          fscanf(dat, "%lf", &Q[i*n_actions*n_theta + j*n_theta + k] );
        }
      }
    }
  }
  else { printf("error loading file\n"); }

  fclose(dat);
}

int check(int s_alpha, int action, int s_theta){

  #ifdef DEBUGSARSA
  //printf("s_alpha %d\n", s_alpha);
  #endif

  int c = 0;
  int matrix_index = s_alpha*n_actions*n_theta + action*n_theta + s_theta;
  if ( (matrix_index >= n_alphas*n_actions*n_theta || matrix_index < 0)){
    printf("s_alpha*n_actions*n_theta + action*n_theta + s_theta = %d\n", matrix_index);
    c = 1;
  }
  return c;
}

void initialize_Q(double *Q, double max_Q_value){

  for (int i=0; i<n_alphas; i++){
    for (int j=0; j<n_actions; j++){
      for (int k=0; k<n_theta; k++){
        if ( (i == 0 && j == 0) || (i==n_alphas-1 && j==n_actions-1)){
          Q[i*n_actions*n_theta + j*n_theta + k] = 0;
        }
        else Q[i*n_actions*n_theta + j*n_theta + k] = max_Q_value;
      }
    }
  }

}

void initialize_Q_counter(int *Q, int value){

  for (int i=0; i<n_alphas; i++){
    for (int j=0; j<n_actions; j++){
      for (int k=0; k<n_theta; k++){
        if ( (i == 0 && j == 0) || (i==n_alphas-1 && j==n_actions-1)){
          Q[i*n_actions*n_theta + j*n_theta + k] = 0;
        }
        else Q[i*n_actions*n_theta + j*n_theta + k] = value;
      }
    }
  }

}

int find_theta_state(double theta){

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

int select_alpha_action(double epsilon, double * Q, int s_alpha, int s_theta, int episode){

  int a1;

  #ifdef DEBUGSARSA
  printf("\n***select action from status = %d\n", s_alpha);
  #endif

  double probab = (double)rand()/RAND_MAX; //rand() % 100;

  if (probab >= epsilon){ // new action as a random number
    #ifdef DEBUGSARSA
        printf("   random\n");
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
  }

  else {

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
      printf("   Q[(s_alpha+%d-1)*n_theta + s_theta] = %f\n", i, Q[(s_alpha+i-1)*n_theta + s_theta]);
      #endif
      if (Q[s_alpha*n_actions*n_theta + i*n_theta + s_theta] > Q[s_alpha*n_actions*n_theta + a1*n_theta + s_theta])  { 
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

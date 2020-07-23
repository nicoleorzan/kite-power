#include <stdio.h>
#include <math.h>

#ifndef __sarsa__
#define __sarsa__

#define PI 3.1415926535897932384626433

#define PENALTY -300.0

#define n_actions 3 // 0 diminuisco alpha, 1 rimango, 2 aumento
#define n_angles 30

#define decision_time 1000
#define max_steps 2000000

#define s_alpha0 10

#ifndef n_alphas
#define n_alphas 15
#endif

double epsilon = 0.1;

double angles[n_angles+1] = {-4, -3.5, -3.2, -3.1, -3.0, -2.9, -2.7, -2.4, -2.1, -1.8, -1.5, -1.2, -0.9, -0.6, \
-0.3, 0., 0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.1, 2.4, 2.7, 2.9, 3.0, 3.1, 3.2, 3.5, 4.}; // impose value n_angles to 30

void load_matrix(double * Q, char * name){

  FILE *dat = fopen(name, "r"); // opens file for reading //"Qmatrix.dat"
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
  else { printf("error loading file\n"); }

  fclose(dat);
}

int check(int s_alpha, int action, int s_vrel_mod){

    int idx = s_vrel_mod*n_alphas*n_actions + s_alpha*n_actions + action;
    int c = 0;
    if ( (idx >= n_alphas*n_actions*n_angles || idx < 0)){
      c = 1;
    }
    return c;
}

int find_state_angle(double v){

  int slice = 0;

  for (int i=0; i<n_angles; i++){
      if (v >= angles[i] && v < angles[i+1]){
          slice = i;
      }
  }

  return slice;
}

void print_matrix(double * Q){

  for (int kk=0; kk<n_angles; kk++){
    printf("angle between %f and %f\n", angles[kk], angles[kk+1]);
    for (int i=0; i<n_alphas; i++){
      printf("alph idx=%d ", i);
      for (int j=0; j<n_actions; j++){
        printf("%f ", Q[kk*n_alphas*n_actions + i*n_actions + j] );
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

int select_greedy_action(double * Q, int s_alpha, int s_vrel_mod){

  int a1;

  /* fix starting action for the search: 0 in any case, and 1 only if the 
    starting state is the beginning of the array */

  int min_value = 0;
  a1 = 0;
  if (s_alpha == 0) { 
    a1 = 1; 
    min_value = 1;
  }

  /* fix ending state: 3 (n_actions), NOT INCLUDED, in any state; 2 (n_actions -1)
    only if starting state is the last one in the array */

  int max_value = n_actions;
  if (s_alpha == n_alphas - 1) { 
    max_value = n_actions - 1; 
  }

  /* LOOP TO SEARCH NEXT ACTION */ 

  for (int i=min_value+1; i<max_value; i++){
    // io sono in Q[s_alpha*n_actions] 
    if (Q[s_vrel_mod*n_alphas*n_actions + s_alpha*n_actions + i] > Q[s_vrel_mod*n_alphas*n_actions + s_alpha*n_actions + a1]) { 
      a1 = i;
    }
  }

  return a1;
}



#endif
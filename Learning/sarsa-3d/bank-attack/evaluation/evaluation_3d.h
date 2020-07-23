#include <stdio.h>
#include <math.h>

#ifndef __sarsa__
#define __sarsa__

#define PI 3.1415926535897932384626433

#define PENALTY -300.0

#define n_actions 3 // 0 diminuisco alpha, 1 rimango, 2 aumento
#define n_bank 13 //7

#define decision_time 1000

#define Gamma 0.9999999999
#define max_steps 10000000

#define s_alpha0 14

#ifndef n_alphas
#define n_alphas 15
#endif

// n bank 13
double bank[n_bank] = {-0.174532925, -0.13962634016, -0.10471975512, -0.06981317008, -0.03490658504, -0.01745329252, 0., 0.01745329252, 0.03490658504, 0.06981317008, 0.10471975512, 0.13962634016, 0.174532925}; //-10, -8, -4, -6, -2, -1, 0, 1, 2, 4, 6, 8, 10
#define s_bank0 6

double epsilon = 0.0001;


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

void select_action_greedy(double * Q, int s_alpha, int *a_alpha, int s_bank, int *a_bank){

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
      }

    }
    }

}


#endif
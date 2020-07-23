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

double angles[n_angles+1] = {-4, -3.5, -3.2, -3.1, -3.0, -2.9, -2.7, -2.4, -2.1, -1.8, -1.5, -1.2, -0.9, -0.6, \
-0.3, 0., 0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.1, 2.4, 2.7, 2.9, 3.0, 3.1, 3.2, 3.5, 4.}; // impose value n_angles to 29

void load_policy(double * pi, char * name){

  FILE *dat = fopen(name, "r"); // opens file for reading
  if (dat)
  {
    for (int kk = 0; kk < n_angles; kk++){
      for (int i = 0; i < n_alphas; i++){
        for (int j = 0; j < n_actions; j++){
          fscanf(dat, "%lf ", &pi[kk*n_alphas*n_actions + i*n_actions + j]);
        }
      }
    }
  }
  else { printf("Nic: error loading file\n"); }

  fclose(dat);
}

void print_policy(double * pi){

    for (int kk = 0; kk < n_angles; kk++){
      printf("angleidx=%d\n", kk);
      for (int i = 0; i < n_alphas; i++){
        for (int j = 0; j < n_actions; j++){
          printf("%f ", pi[kk*n_alphas*n_actions + i*n_actions + j]);
        }
        printf("\n");
      }
      printf("\n");
    }

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

int select_action(double * pi, int s_alpha, int s_vrel_mod){

  int a1;
  double a, b, c;
  a = pi[s_vrel_mod*n_alphas*n_actions + s_alpha*n_actions + 0];
  b = pi[s_vrel_mod*n_alphas*n_actions + s_alpha*n_actions + 1];
  c = pi[s_vrel_mod*n_alphas*n_actions + s_alpha*n_actions + 2];

  double probab = (double)rand()/RAND_MAX;

  if (probab >= 0 && probab < a){
    a1 = 0;
  }
  else if (probab >= a && probab < a+b){
    a1 = 1;
  }
  else if (probab >= a+b && probab < a+b+c){
    a1 = 2;
  }

  return a1;
}



#endif
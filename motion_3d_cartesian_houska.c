#include "Dynamics/dynamics_3d_cartesian_houska.h"
#include "Dynamics/winds.h"
#include <time.h>
#include <string.h>
#include <stdbool.h>

#define _theta0 PI/4.
#define _phi0 0.
#define _dtheta0 0.
#define _dphi0 0.
#define dim 3
#define mu 0.261799//=15 //0.0872665=5          RADIANS!!!

// ============== FILE INPUT: ATTACK ANGLE AND WIND X, Y, Z, THETA0 AND PHI0 ============

int main(int argc, char *argv[]){

    // ========================= READING INPUT VARIABLES ==========================

    if (argc <= 4){
        printf("Missing inputs!!! Need ATTACK ANGLE INDEX, WIND X, WIND Y, WIND Z, THETA0 AND PHI0\n");
        return 0;
    }

    int alpha_index = atoi( *(argv + 1) );

    double W[3];
    W[0] = atof( *(argv + 2) );
    W[1] = atof( *(argv + 3) );
    W[2] = atof( *(argv + 4) );

    double theta = _theta0;
    double phi = _phi0;
    double dtheta = _dtheta0;
    double dphi = _dphi0;

    int uno=0, due=0, tre=0;
    double sp=0, pv=0;

    if (alpha_index >= n_alphas){
        printf("Alpha index too big!!!\n");
        return 0;
    }

    // ========================= OUTPUT FILES ==========================

    FILE *trajectory, *debug;
    trajectory = fopen("out.txt", "w+");
    debug = fopen("hdebug3d.csv", "w+");

    //fprintf(trajectory, "t,x_kite,y_kite,z_kite,x_block,y_block,z_block,theta,vtheta,windx,windy,wind_z,v_blockx,v_blocky,Tension\n");
    //fprintf(debug, "i,Alpha,mu,theta,Windx,Windy,Windz,Vkx,Vky,Vkz,Lift,Liftx,Lifty,Liftz,Drag,Tension,F_attrito,sector,uno,due,tre,1-sp,angle,prod_vect,t2[0],t2[1],t2[2]\n");

    // ============================ VARIABLES DEFINITION ============================

    // kite motion vectors from fixed origin (x, z)
    double *rk = (double*) malloc(dim * sizeof(double)); 
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

    double *t22= (double*) malloc(3 * sizeof(double)); 

    double r_diff_modulo;
    double v_diff_modulo;

    double lift=0, drag=0;
    double T = 0;
    double F_attr;
    double F_vinc;
    double theta_star;
    int stability = 0;
    int decollato = 0;

    int sector = 0;

    double l0, l1, l2;
    double d0, d1, d2;

    variables_initialization(rk, vk, ak, theta, phi, dtheta, dphi, r_block, v_block, a_block, r_diff, v_diff, a_diff);
    
    //streamfunction3d(rk, W);

    int t = 0;   

    // ============================ KITE FLYING LOOP, STOPS WHEN FALL OCCURS ============================
    
    for (int i=0; i<STEPS; i++){

        integration_trajectory(rk, vk, ak, r_block, v_block, a_block, r_diff, v_diff, a_diff, \
                            &theta, &phi, alpha_index, mu, W, &lift, &drag, &T, &F_attr, i, &sector, \
                            &l0, &l1, &l2, &d0, &d1, &d2);//, &uno, &due, &tre, &sp, &pv, t22);

        //printf("\ni=%d, L=%f, D=%f, T=%f, F_attr=%f, sector=%d\n", i, lift, drag, T, fabs(F_attr), sector);
        //printf("Lx=%f, Ly=%f, Lz=%f\n", l0, l1, l2);

        r_diff_modulo = sqrt(r_diff[0]*r_diff[0] + r_diff[1]*r_diff[1] + r_diff[2]*r_diff[2]);
        F_vinc = m_block*g - T*sin(theta);

        // moving the kite to put it again at distance R with the block

        rk[0] = r_block[0] + (rk[0] - r_block[0])/fabs(r_diff_modulo)*R;
        rk[1] = r_block[1] + (rk[1] - r_block[1])/fabs(r_diff_modulo)*R;
        rk[2] = r_block[2] + (rk[2] - r_block[2])/fabs(r_diff_modulo)*R;

        //streamfunction3d(rk, W);

        /*if (i%50 == 0){
            fprintf(debug, "%d,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%d,%d,%d,%d,%f,%f,%f,%f,%f,%f\n", \
                    i, alphas[alpha_index], mu, theta, W[0], W[1], W[2], vk[0], vk[1], vk[2], lift, l0, l1, l2, drag, T, \
                    fabs(F_attr), sector, uno, due, tre, 1-fabs(sp), acos(sp), pv, t22[0], t22[1], t22[2]);
        }*/

        if (i%PRINTSTEP == 0 || rk[2] <= 0.){
            //fprintf(trajectory, "%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", \
            //        t, rk[0], rk[1], rk[2], r_block[0], r_block[1],r_block[2], theta, dtheta, \
                    W[0], W[1], W[2], v_block[0], v_block[1],T);
            if (rk[2] <=0. ){
                //printf("Kite Fall, steps %d, z<0, break\n", i);
                break;
            }
        }

        t += 1;

        if (F_vinc < 0) {
            decollato = 1;
        }

    }

    v_diff_modulo = sqrt(v_diff[0]*v_diff[0] + v_diff[1]*v_diff[1] + v_diff[2]*v_diff[2]);

    dtheta = -1/(sqrt(1-(r_diff[2]/r_diff_modulo)*(r_diff[2]/r_diff_modulo)))*(v_diff[2]*r_diff_modulo-r_diff[2]*v_diff_modulo)/(r_diff_modulo*r_diff_modulo);

    dphi = 1/(1+(r_diff[1]/r_diff[0])*(r_diff[1]/r_diff[0]))*(v_diff[1]*r_diff[0] - r_diff[0]*v_diff[1])/(r_diff[0]*r_diff[0]);
    
    if (rk[2] <= 0){
        rk[2] = 0;
        v_block[0] = 0;
        v_block[1] = 0;
        vk[0] = 0;
        vk[1] = 0;
        vk[2] = 0;
        theta = PI/2;
        dtheta = 0;
    }

    //printf("\niter, tot time, m_block, alpha, mu, theta0, theta_fin, v_block_fin_x, v_block_fin_y, Wind_x, Wind_y, Wind_z");
    //printf(" Vkitex, Vkitey, Vkitez, vrelkite_x, vrelkite_y, vrelkite_z, F_vinc, Tension, Lift, l0, l1, l2, Drag, d0, d1, d2\n");
    
    printf("%d, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f,  %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f\n", \
    t, t*h, m_block, alphas[alpha_index], mu, _theta0, theta, dtheta, v_block[0], v_block[1], W[0], W[1], W[2], \
    vk[0], vk[1], vk[2], vk[0] - W[0], vk[1] - W[1], vk[2] - W[2], F_vinc, T, lift, l0, l1, l2, drag, d0, d1, d2);

    free(rk);
    free(vk);
    free(ak);

    free(r_block);
    free(v_block);
    free(a_block);

    free(r_diff);
    free(v_diff);
    free(a_diff);

    fclose(trajectory);
    fclose(debug);

    //remove("a.out");

    return 0;

    }

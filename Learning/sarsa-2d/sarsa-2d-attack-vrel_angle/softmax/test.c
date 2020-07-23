#include "../../../../Dynamics/dynamics_2d_cartesian.h"
#include "../../../../Dynamics/winds.h"
#include "softmax.h"

double reward_dt = h * decision_time;
#define dim 2

int main(int argc, char *argv[])
{

    FILE *out, *infos;
    out = fopen("cout.txt", "w");
    infos = fopen("infos.txt", "w");
    fprintf(out, "step,x_kite,z_kite,x_block,z_block,theta,windx,windy,v_block,vrelx,vrely,v_rel_mod,vrel_angle,s_vrel_angle,return,alpha,action,vrel_angle_idx,reward\n");
    fprintf(infos, "relative velocities partition\n");

    for (int i = 0; i < n_angles + 1; i++)
    {
        fprintf(infos, "%f, ", angles[i]);
    }
    fprintf(infos, "\n");
    fclose(infos);

    // ======== DYNAMICS VARIABLES =======

    // vettori moto kite dall'origine fissa (x, z)
    double *rk = (double *)malloc(dim * sizeof(double));
    double *vk = (double *)malloc(dim * sizeof(double));
    double *ak = (double *)malloc(dim * sizeof(double));

    // vettori moto block dall'origine fissa (x, z)
    double *r_block = (double *)malloc(dim * sizeof(double));
    double *v_block = (double *)malloc(dim * sizeof(double));
    double *a_block = (double *)malloc(dim * sizeof(double));

    double *r_diff = (double *)malloc(dim * sizeof(double));
    double *v_diff = (double *)malloc(dim * sizeof(double));
    double *a_diff = (double *)malloc(dim * sizeof(double));

    double r_diff_modulo;

    double theta;
    double vtheta;

    double W[dim] = {0, 0};
    int et_val = 0;

    //double lr = 1.0; // adaptive learning rate
    double vrel_x, vrel_z;
    double vrel_angle;

    double T = 0;
    double F_attr = 0;
    double lift = 0, drag = 0;
    int sector = 0;
    double l0, l1;
    double d0, d1;

    // ======== LEARNING VARIABLES =======

    double *pi = (double *)malloc(n_alphas * n_actions * n_angles * sizeof(double));

    double reward = 0.;
    double tot_reward = 0.;

    // SARSA STATE
    int s_alpha, s_alpha1;
    int s_vrel_angle, s_vrel_angle1;

    // SARSA ACTION
    int a_alpha = 0, a_alpha1 = 0;

    int it = 0;

    load_policy(pi, "softmaxpolicy.txt");
    print_policy(pi);

    // ======================= INITIALIZATION ==========================

    for (int episode=0; episode<20; episode ++){

        printf("Policy test\n");

        theta = theta0;
        vtheta = vtheta0;

        variables_initialization(rk, vk, ak, theta, vtheta, r_block, v_block, a_block, r_diff, v_diff, a_diff);

        streamfunction2d_hard(rk, W);

        s_alpha = s_alpha0;

        vrel_x = vk[0] - W[0];
        vrel_z = vk[1] - W[1];

        vrel_angle = atan2(vrel_z, vrel_x);
        s_vrel_angle = find_state_angle(vrel_angle);

        //printf("vrelx=%f, vrelz=%f, angle=%f, angle state=%d\n", vrel_x, vrel_z, vrel_angle, s_vrel_angle);

        a_alpha = select_action(pi, s_alpha, s_vrel_angle);

        it = 0;
        reward = 0;
        tot_reward = 0.;

        printf("theta0 = %f, vel_theta0 = %f\n", theta, vtheta);
        printf("initial s_alpha=%d\n", s_alpha);
        printf("W[0]=%f, W[1]=%f\n", W[0], W[1]);

        integration_trajectory(rk, vk, ak, r_block, v_block, a_block, r_diff, v_diff, a_diff,
                            &theta, &vtheta, s_alpha, W, &lift, &drag, &T, &F_attr, it, &sector, &l0, &l1, &d0, &d1, &et_val);

        streamfunction2d_hard(rk, W);

        //printf("s_alpha=%d, a_alpha=%d, s_vrel=%d\n", s_alpha, a_alpha, s_vrel_angle);

        while (rk[1] > 0)
        {

            if (it % decision_time == 0){

                //printf("s_alpha=%d, a_alpha=%d, s_vrel=%d\n", s_alpha, a_alpha, s_vrel_angle);

                fprintf(out, "%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%d,%f,%f,%d,%d,%f\n", it, rk[0], rk[1],
                        r_block[0], r_block[1], theta, W[0], W[1], v_block[0], vrel_x, vrel_z,
                        sqrt(vrel_x * vrel_x + vrel_z * vrel_z), vrel_angle, s_vrel_angle, tot_reward, alphas[s_alpha], a_alpha, s_vrel_angle, reward);
            }

            if (it > max_steps) {
                printf("MAX STEPS, %d, exiting\n", max_steps);
                printf("final alpha=%d\n", s_alpha);
                printf("return=%f, space percurred=%f\n\n", tot_reward, r_block[0]);
                break;
            }

            integration_trajectory(rk, vk, ak, r_block, v_block, a_block, r_diff, v_diff, a_diff,
                                &theta, &vtheta, s_alpha, W, &lift, &drag, &T, &F_attr, it, &sector, &l0, &l1, &d0, &d1, &et_val);

            r_diff_modulo = sqrt(r_diff[0] * r_diff[0] + r_diff[1] * r_diff[1]);

            rk[0] = r_block[0] + (rk[0] - r_block[0]) / fabs(r_diff_modulo) * R;
            rk[1] = r_block[1] + (rk[1] - r_block[1]) / fabs(r_diff_modulo) * R;

            streamfunction2d_hard(rk, W);

            reward = fabs(v_block[0]) * reward_dt;

            if (it % decision_time == 0){
                tot_reward += reward;
            }

            if (rk[1] <= 0.){
                printf("Kite fallen: z<0, steps=%d, break\n", it);
                printf("final alpha=%d\n\n", s_alpha);
                
                fprintf(out, "%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%d,%f,%f,%d,%d,%f\n", it, rk[0], rk[1],
                        r_block[0], r_block[1], theta, W[0], W[1], v_block[0], vrel_x, vrel_z,
                        sqrt(vrel_x * vrel_x + vrel_z * vrel_z), vrel_angle, s_vrel_angle, tot_reward, alphas[s_alpha], a_alpha, s_vrel_angle, reward);
                break;
            }

            if (it % decision_time == 0){

                // FIND STATE S1 USING ACTION A

                s_alpha1 = update_state(s_alpha, a_alpha);

                vrel_x = vk[0] - W[0];
                vrel_z = vk[1] - W[1];
                vrel_angle = atan2(vrel_z, vrel_x);
                s_vrel_angle1 = find_state_angle(vrel_angle);

                // SEARCH FOR NEXT ACTION A1

                a_alpha1 = select_action(pi, s_alpha1, s_vrel_angle1);

                // MOVE ON: S = S1, A = A1
                s_alpha = s_alpha1;
                s_vrel_angle = s_vrel_angle1;

                a_alpha = a_alpha1;
            }

            it += 1;
        }

        episode +=1;

    }

    free(rk);
    free(vk);
    free(ak);

    free(r_diff);
    free(v_diff);
    free(a_diff);

    free(r_block);
    free(v_block);
    free(a_block);

    free(pi);

    fclose(out);

    return 0;
}

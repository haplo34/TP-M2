/* Algorithme velocity-verlet avec un potentiel de Lenard-Jones
 * pour un système de N particules */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#define EPSILON 1
#define SIGMA 1


/* Defining factorial function to determine the number
 * of interactions. */
long factorial(int n)
{
    if (n==0)
        return 1;
    else
        return n*factorial(n-1);
}

    
/* Lennard-Jones potential (force) */
void lennard_jones_f(double x1, double x2, double y1, double y2, double *fx, double *fy)
{
    double dx = x2 - x1;
    double dy = y2 - y1;
    double r2 = dx * dx + dy * dy;
    
    *fx += EPSILON * (48 * dx * pow(SIGMA, 12) / pow(r2, 7)
    - 24 * dx * pow(SIGMA, 6) / pow  (r2, 4));
    
    *fy += EPSILON * (48 * dy * pow(SIGMA, 12) / pow(r2, 7)
    - 24 * dy * pow(SIGMA, 6) / pow  (r2, 4));
}


/* Velocity-Verlet algoritm 
 * p[i][0] = x pour la particule i, p[i][1] = y
 * p[i][2] = v_x, p[i][3] = v_y
 */
void velocity_verlet(double p[][4], int Np, double m, double dt)
{
    for (int i=0; i<Np; i++)
    {
        double fx_0 = 0; double fy_0 = 0;
        double fx_1 = 0; double fy_1 = 0;
        
        for (int j=0; j<Np; j++)
        {
            if (i != j)
            {
                lennard_jones_f(
                    p[i][0], p[j][0], p[i][1], p[j][1], &fx_0, &fy_0);
            }
        }
        p[i][0] += p[i][2] * dt + 1/(2*m) * fx_0 * dt * dt;
        p[i][1] += p[i][3] * dt + 1/(2*m) * fy_0 * dt * dt;
        
        for (int j=0; j<Np; j++)
        {
            if (i != j)
            {
                lennard_jones_f(
                    p[i][0], p[j][0], p[i][1], p[j][1], &fx_1, &fy_1);
            }
        }
        p[i][2] += 1/(2*m) * (fx_1 + fx_0) * dt;
        p[i][3] += 1/(2*m) * (fy_1 + fy_0) * dt;
    }
}


int main(int argc, char *argv[])
{
    // 1st argument is number of particules
    int Np = atoi(argv[1]);
    // 2nd argument is mass of particules
    double m = strtod(argv[2], NULL);
    // 3rd argument is step of computation
    double dt = strtod(argv[3], NULL);
    // 4th argument is stop condition
    double t_max = strtod(argv[4], NULL);
    // 2D array declaration containing the particules
    double p[][4] = 
    {
        {1, 5, 0, 0},
        {50, 25, 0, 0},
        {100, 15, 0, 0},
        {35, 200, 0, 0}
    };
    // Number of interactions
    int Ni = factorial(Np - 1);
    
    // Velocity-Verlet call
    double t = 0;
    
    while (t < t_max)
    {
        t += dt;
        velocity_verlet(p, Np, m, dt);
    }
    
    for (int i=0; i<Np; i++)
    {
        for (int j=0; j<4; j++)
        {
            printf("%f  ", p[i][j]);
        }
        printf("\n");
    }
    return 0;
}

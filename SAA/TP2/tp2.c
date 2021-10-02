#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define EPSILON 1
#define SIGMA 1
const double Rc=2.5;  /* global constant */


void interactions(double x1, double x2, double y1, double y2,
                     double *u12, double *fx, double *fy)
{
    double dx = x2 - x1;
    double dy = y2 - y1;
    double r2 = dx * dx + dy * dy;
    double r6 = r2 * r2 * r2;

    *u12 += 4 * EPSILON * (
        pow(SIGMA, 12) / (r6 * r6) - pow(SIGMA, 6) / r6);
    
    *fx += EPSILON * (48 * dx * pow(SIGMA, 12) / pow(r2, 7)
        - 24 * dx * pow(SIGMA, 6) / pow  (r2, 4));
    
    *fy += EPSILON * (48 * dy * pow(SIGMA, 12) / pow(r2, 7)
        - 24 * dy * pow(SIGMA, 6) / pow  (r2, 4));
}


int main()
{
    double u12, fx, fy;
    double p[][4] = 
    {
        {1, 5, 0, 0},
        {50, 25, 0, 0},
    };
    interactions(p[0][0], p[1][0], p[0][1], p[1][1], &u12, &fx, &fy);
    return 0;
}
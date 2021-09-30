#include <stdio.h>

#define EPSILON 1
#define SIGMA 1
const double Rc=2.5  /* global constant */




void interactions(double x1, double y1, double x2, double y2,
                  double &u12, double &f12x, double &f12y)
{
    double dx = x2 - x1;
    double dy = y2 - y1;
    double r2 = pow(dx, 2) + pow(dy, 2);
    
    u12 = 4 * EPSION * ((SIGMA**2 / r2)**6 - (SIGMA**2 / r2)**3)
    
    f12x = EPSILON * (48 * dx * pow(SIGMA, 12) / pow(r2, 7)
    - 24 * dx * pow(SIGMA, 6) / pow  (r2, 4));
    
    f12y = EPSILON * (48 * dy * pow(SIGMA, 12) / pow(r2, 7)
    - 24 * dy * pow(SIGMA, 6) / pow  (r2, 4));
}

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double deriv_f(double t, double var) {
    return exp(var - t)*(1.0/cos(var))*(1.0 + t*t);
}

double doubl_div(double a, double b) {
    return a/b;
}

int main() {
    int i, j=0, m=7, n=152;
    double *k, RKF_4, RKF_5, y_n=0.0, t_n=0.0, step=0.25, err, TOL=5.0*1.0e-10;
    if ((k = (double *) malloc(m*sizeof(double))) == NULL) {
        fprintf(stderr,"malloc failed\n");
        exit(1);
    }
    fprintf(stderr, "i   t_n                y_n                step\n");
    fprintf(stderr, "%d,  %.12f,    %.12f,    %.12f\n", i, t_n, y_n, step);
    t_n += step;
    for (i=1; i<n; i++) {
        k[0] = deriv_f(t_n, y_n);
        k[1] = deriv_f(t_n+doubl_div(1.0,5.0)*step, y_n+doubl_div(1.0,5.0)*step*k[0]);
        k[2] = deriv_f(t_n+doubl_div(3.0,10.0)*step, y_n+doubl_div(3.0,40.0)*step*k[0]+doubl_div(9.0,40.0)*step*k[1]);
        k[3] = deriv_f(t_n+doubl_div(4.0,5.0)*step, y_n+doubl_div(44.0,45.0)*step*k[0]-doubl_div(56.0,15.0)*step*k[1]+doubl_div(32.0,9.0)*step*k[2]);
        k[4] = deriv_f(t_n+doubl_div(8.0,9.0)*step, y_n+doubl_div(19372.0,6561.0)*step*k[0]-doubl_div(25360.0,2187.0)*step*k[1]+doubl_div(64448.0,6561.0)*step*k[2]-doubl_div(212.0,729.0)*step*k[3]);
        k[5] = deriv_f(t_n+step, y_n+doubl_div(9017.0,3168.0)*step*k[0]-doubl_div(355.0,33.0)*step*k[1]+doubl_div(46732.0,5247.0)*step*k[2]+doubl_div(49.0,176.0)*step*k[3]-doubl_div(5103.0,18656.0)*step*k[4]);
        k[6] = deriv_f(t_n+step, y_n+doubl_div(35.0,384.0)*step*k[0]+doubl_div(500.0,1113.0)*step*k[2]+doubl_div(125.0,192.0)*step*k[3]-doubl_div(2187.0,6784.0)*step*k[4]+doubl_div(11.0,84.0)*step*k[5]);
        RKF_4 = y_n + step*(doubl_div(5179.0,57600.0)*k[0]+doubl_div(7571.0,16695.0)*k[2]+doubl_div(393.0,640.0)*k[3]-doubl_div(92097.0,339200.0)*k[4]+doubl_div(187.0,2100.0)*k[5]+doubl_div(1.0,40.0)*k[6]);
        RKF_5 = y_n + step*(doubl_div(35.0,384.0)*k[0]+doubl_div(500.0,1113.0)*k[2]+doubl_div(125.0,192.0)*k[3]-doubl_div(2187.0,6784.0)*k[4]+doubl_div(11.0,84.0)*k[5]);
        err = fabs(RKF_4 - RKF_5);
        //fprintf(stderr, "%d, %f\n" , i, err);
        if (err >= doubl_div(1.0,4.0)*step*TOL && err <= step*TOL) {
            y_n = RKF_5;
            fprintf(stderr, "%d,  %.12f,    %.12f,    %.12f\n", i, t_n, y_n, step);
            t_n += step;
        }
        else if (err > step*TOL) {
            step /= (double)2.0;
            i -= 1;
            t_n -= step;
        }
        else if (err < doubl_div(1.0,4.0)*step*TOL) {
            y_n = RKF_5;
            fprintf(stderr, "%d,  %.12f,    %.12f,    %.12f\n", i, t_n, y_n, step);
            step *= (double)2.0;
            t_n += step;
        }
    }
    free(k);
}
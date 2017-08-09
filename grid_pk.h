#ifndef GRID_PK_H
#define GRID_PK_H
#include "grid_fft.h"
#include "read_athena_tracers.hpp"

//double *grid_tsc(double *x, double *y, double *z, double *m, int N, FFTW_Grid_Info grid_info);
double *grid_tsc(vector<tracer>t, FFTW_Grid_Info grid_info);
double *grid_cic(double *x, double *y, double *z, double *m, int N, FFTW_Grid_Info grid_info);
double *grid_ngp(double *x, double *y, double *z, double *m, int N, FFTW_Grid_Info grid_info);
double *grid_dfk(int N, double *u, FFTW_Grid_Info grid_info, MPI_Comm world);
double w_p(int p, double kx, double ky, double kz, FFTW_Grid_Info grid_info);
#endif //GRID_PK_H

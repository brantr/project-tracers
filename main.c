#include <stdio.h>
#include <math.h>
#include <vector>
#include "rng.h"
#include "grid_fft.h"
#include "grid_pk.h"
#include "read_athena_tracers.hpp"
struct Particle
{
  double *m;
  double *x;
  double *y;
  double *z;
};


int main(int argc, char **argv)
{
  int N = 100000; //number of objects

  struct Particle p;

  FILE *fp;

  FFTW_Grid_Info grid_info;
  char fname[200];
  char filename[200];
  char fdir[200];
  char fsuff[200];


  int myid;
  int numprocs;
  int isnap;

  double *u;
  double *dfk;
  double A     = 1.0;
  double sigma = 0.1;
  int ns = 10;
  int nx = 64;
  int ny = 64;
  int nz = 64;
  int nfiles;
  int pmode = 3;

  AthenaHeaderH *h;

  //check number of arguments
  if(!((argc==3)||(argc==5)||(argc==6)||(argc==7)||(argc==8)))
  {
    if(myid==0)
    {
      printf("./athena4_binary.tracer_projection fdir filename [isnap nfiles] {suff} ngrid mode\n");
      fflush(stdout);
    }
    MPI_Finalize();
    exit(-1);
  }

  //print the data directory and filename
  sprintf(fdir,"%s",argv[1]);
  sprintf(filename,"%s",argv[2]);

  //check to see if the binary files
  //are in a monolithic binary
  if(argc>=5)
  {
    isnap    = atoi(argv[3]);
    nfiles = atoi(argv[4]);
  }else{
    nfiles = 1;
  }
  if(argc>=6)
    sprintf(fsuff,"%s",argv[5]);
  if(argc>=7)
  {
    nx = atoi(argv[6]);
    ny = nx;
    nz = nx;
  }
  if(argc==8)
    pmode = atoi(argv[7]);

  printf("Directory = %s.\n",fdir);
  printf("Filename  = %s.\n",filename);
  printf("isnap     = %d.\n",isnap);
  printf("nfiles    = %d.\n",nfiles);
  if(argc>=6)
    printf("suffix    = %s.\n",fsuff);


  grid_info.nx = nx;
  grid_info.ny = ny;
  grid_info.nz = nz;
  grid_info.ndim = 3;
  grid_info.BoxSize = 1.0;


  char fname_tracers[200];
  long nt;
  int count=0;
  int isub = 0;
  float xmin;
  float xmax;
  FILE *fpx;

  vector<tracer> t;
  vector<tracer> tin;


  printf("ngrid = %d\n",nx);
  printf("pmode = %d\n",pmode);

  MPI_Init(&argc,&argv);
  MPI_Comm world = MPI_COMM_WORLD;
  MPI_Comm_rank(world,&myid);
  MPI_Comm_size(world,&numprocs);



  //set the grid size
  initialize_mpi_local_sizes(&grid_info, world);

  h = NULL;

  //Loop over the number of files
  //for(isub=myid;isub<nfiles;isub+=numprocs)
  //for(isub=0;isub<1;isub++)
  for(isub=myid;isub<nfiles;isub+=numprocs)
  {
    if(isub==0)
    {
      //create a new filename
      sprintf(fname_tracers,"%s/id%d/%s.%04d.tra",fdir,isub,filename,isnap);
    }else{
      //create a new filename
      sprintf(fname_tracers,"%s/id%d/%s-id%d.%04d.tra",fdir,isub,filename,isub,isnap);
    }

    //read in the tracers from the file
    if(h)
      free(h);
    h = read_athena_tracers(fname_tracers,&tin);



    t.insert(t.end(), tin.begin(), tin.end());
    //print information about the tracers
    printf("myid %2d isub %3d tin.size() %10ld t.size() %10ld\n",myid,isub,tin.size(),t.size());
    vector<tracer>().swap(tin);

  }
  printf("total number of tracers = %ld\n",t.size());


  nx = h->nx;
  ny = h->ny;
  nz = h->nz;

  ny = 512;
  nz = 512;
  double *xout = (double *) calloc(ny*nz,sizeof(double *));

  int ii, jj;
  double x1, x2;
  for(long i=0;i<t.size();i++)
  {
    x1 = t[i].x[1];
    x2 = t[i].x[2];

    if(x1<0)
      x1 += 1.;
    if(x1>=1.)
      x1 -= 1.0;
    if(x2<0)
      x2 += 1.0;
    if(x2>=1.0)
      x2 -= 1.0;
    ii = (int) ((double) ny*x1);
    jj = (int) ((double) nz*x2);

    if(i<10)
      printf("x1 %e x2 %e ii %d jj %d\n",x1,x2,ii,jj);
    xout[ii*nz + jj] += 1.0;

  }

  char fxout[200];
  sprintf(fxout,"tracers.%04d.dat",isnap);
  fp = fopen(fxout,"w");
  fwrite(&ny,sizeof(int),1,fp);
  fwrite(&nz,sizeof(int),1,fp);
  fwrite(xout,sizeof(double),ny*nz,fp);
  fclose(fp);
  for(long i=0;i<t.size();i++)
  {
    x1 = t[i].x[1];
    x2 = t[i].x[2];

    if(x1<0)
      x1 += 1.;
    if(x1>=1.)
      x1 -= 1.0;
    if(x2<0)
      x2 += 1.0;
    if(x2>=1.0)
      x2 -= 1.0;
    ii = (int) ((double) ny*x1);
    jj = (int) ((double) nz*x2);

    if(i<10)
      printf("x1 %e x2 %e ii %d jj %d\n",x1,x2,ii,jj);
    xout[ii*nz + jj] += t[i].d;

  }
  sprintf(fxout,"tracers.%04d.density.dat",isnap);
  fp = fopen(fxout,"w");
  fwrite(&ny,sizeof(int),1,fp);
  fwrite(&nz,sizeof(int),1,fp);
  fwrite(xout,sizeof(double),ny*nz,fp);
  fclose(fp);
  free(xout);
  MPI_Finalize();
  return 0;

  //make sure particles are in 0<x<1
  double xxmin[3] = {1.0e9,1.0e9,1.0e9};
  double xxmax[3] = {-1.0e9,-1.0e9,-1.0e9};
  for(long i=0;i<t.size();i++)
  {
    t[i].m = 1.0;
    for(int k=0;k<3;k++)
    {
      if(t[i].x[k]<0)
        t[i].x[k] += 1.0;
      if(t[i].x[k]>=1.0)
        t[i].x[k] -= 1.0;

      if(t[i].x[k]<xxmin[k])
        xxmin[k] = t[i].x[k];
      if(t[i].x[k]>xxmax[k])
        xxmax[k] = t[i].x[k];
    }
  }
  for(int k=0;k<3;k++)
    printf("xmin %e xmax %e\n",xxmin[k],xxmax[k]);

  printf("nx %ld ny %ld nz %ld BoxSize %f\n",grid_info.nx,grid_info.ny,grid_info.nz,grid_info.BoxSize);

  printf("Performing interpolation on to the grid...\n");

  //at this point, the tracers are all loaded into the array.
  //so we just need to grid and compute as before
  u = grid_tsc(t,grid_info);

  printf("...done!\n");

  //output the ngp grid
  switch(pmode)
  {
    case 1:
      sprintf(fname,"ngp.%s.%04d.%d.dat",filename,isnap,nx);
      break;
    case 2:
      sprintf(fname,"cic.%s.%04d.%d.dat",filename,isnap,nx);
      break;
    default:
      sprintf(fname,"tsc.%s.%04d.%d.dat",filename,isnap,nx);
      break;
  }
  output_fft_grid(fname, u, grid_info, 0, nx, 0, ny, 0, nz, myid, numprocs, world);
  printf("output %s..\n",fname);

  //compute <|df(k)|^2>
  dfk = grid_dfk(1,u,grid_info,world);

  //output the grid
  sprintf(fname,"dfk.%s.%04d.%d.dat",filename,isnap,nx);
  output_fft_grid(fname, dfk, grid_info, 0, nx, 0, ny, 0, nz, myid, numprocs, world);

  //remember numer of particles
  sprintf(fname,"npart.%s.%04d.%d.txt",filename,isnap,nx);
  fp = fopen(fname,"w");
  fprintf(fp,"%ld\n",t.size());
  fclose(fp);

  free(u);
  free(dfk);
  MPI_Finalize();
  exit(0);

/*
  ns = 1;


  initialize_mpi_local_sizes(&grid_info, world);
  for(int id=0;id<ns;id++)
  {

    //create_normal_distribution(&p,N,sigma);
    //create_uniform_distribution(&p,N,A,id);
    read_particle_distribution(&p, &N, fname);
    printf("N = %d\n",N);


    //printf("p[0] %e %e %e\n",p.x[0],p.y[0],p.z[0]);

    //grid particles
    u = grid_ngp(p.x,p.y,p.z,p.m,N,grid_info);
    //u = grid_cic(p.x,p.y,p.z,p.m,N,grid_info);
    //u = grid_tsc(p.x,p.y,p.z,p.m,N,grid_info);

    //u =  grid_make_gaussian_kernel(A, sigma*grid_info.nx, grid_info);

    //output the ngp grid
    sprintf(fname,"ngp.%d.dat",id);
    output_fft_grid(fname, u, grid_info, 0, nx, 0, ny, 0, nz, myid, numprocs, world);
    printf("output...\n");

    //compute <|df(k)|^2>
    dfk = grid_dfk(N,u,grid_info,world);

    //output the grid
    sprintf(fname,"dfk.%d.dat",id);
    output_fft_grid(fname, dfk, grid_info, 0, nx, 0, ny, 0, nz, myid, numprocs, world);

    

    if(id==0)
    {
      FILE *fp;
      sprintf(fname,"particles.%d.dat",id);
      fp = fopen(fname,"w");
      for(int i=0;i<N;i++)
        fprintf(fp,"%e\t%e\t%e\n",p.x[i],p.y[i],p.z[i]);
      fclose(fp);
    }

    free(p.x);
    free(p.y);
    free(p.z);
    free(p.m);

    free(u);
    free(dfk);
  }
  return 0;
*/
}

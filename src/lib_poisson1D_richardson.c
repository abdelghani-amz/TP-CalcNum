/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

void eig_poisson1D(double* eigval, int *la){

  for (size_t ii = 0; ii < *la; ii++) {
    double scal = (1.0 * ii + 1.0) * M_PI_2 * (1.0 / (*la + 1));
    eigval[ii] = 4 * sin(scal) * sin(scal);
  } 
}

double eigmax_poisson1D(int *la){
  double eigmax = sin(*la * M_PI_2 * (1.0 / (*la + 1)));
  return 4 * eigmax * eigmax;
}

double eigmin_poisson1D(int *la){
  double eigmin = sin(M_PI_2 * (1.0 / (*la + 1)));
  return 4 * eigmin * eigmin;
}

double richardson_alpha_opt(int *la){
  return 2 / (eigmax_poisson1D(la) + eigmin_poisson1D(la));
}

void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
  double *tmp = malloc(sizeof(double)*(*la));
  cblas_dcopy(*la,RHS,1,tmp,1);
  cblas_dgbmv(CblasColMajor,CblasNoTrans,*la,*la,*kl,*ku,-1.0,AB,*lab,X,1,1.0,tmp,1);

  double check = cblas_dnrm2(*la,tmp,1);
  *nbite=0;
  resvec[*nbite] = check;

  while(check > (*tol) && *nbite < *maxit)
  {
    cblas_dcopy(*la,RHS,1,tmp,1);
    cblas_dgbmv(CblasColMajor,CblasNoTrans,*la,*la,*kl,*ku,-1.0,AB,*lab,X,1,1.0,tmp,1);
    cblas_daxpy(*la,(*alpha_rich),tmp,1,X,1);
    check = cblas_dnrm2(*la,tmp,1);
    (*nbite)++;
    resvec[*nbite] = check;
  }
  free(tmp);
}

void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
  int ld;
  for(int i = 0; i<(*la); i++)
  {
    ld = i*(*lab);
    MB[ld+*kv] = AB[ld+*kv];
  }
}

void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
  int ld;
  for(int i = 0; i<(*la); i++)
  {
    ld = i*(*lab);
    MB[ld+*kv] = AB[ld+*kv];
    MB[ld+*kv+1] = AB[ld+*kv+1];
  }
}

void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
  double *tmp = malloc(sizeof(double)*(*la));
  cblas_dcopy(*la,RHS,1,tmp,1);
  cblas_dgbmv(CblasColMajor,CblasNoTrans,*la,*la,*kl,*ku,-1.0,AB,*lab,X,1,1.0,tmp,1);

  int info; int NHRS = 1;
  int *ipiv = (int *)calloc(*la, sizeof(int));
  int ku_mb = (*ku)-1;
  dgbtrf_(la,la,kl,&ku_mb,MB,lab,ipiv,&info);


  double check = cblas_dnrm2(*la,tmp,1);
  *nbite=0;
  resvec[*nbite] = check;
  
  while((check > (*tol)) && ((*nbite) < (*maxit)))
  {
    cblas_dcopy(*la,RHS,1,tmp,1);
    cblas_dgbmv(CblasColMajor,CblasNoTrans,*la,*la,*kl,*ku,-1.0,AB,*lab,X,1,1.0,tmp,1);
    dgbtrs_("N",la,kl,&ku_mb,&NHRS,MB,lab,ipiv,tmp,la,&info);
    cblas_daxpy(*la,1.0,tmp,1,X,1);

    check = cblas_dnrm2(*la,tmp,1);
    (*nbite)++;
    resvec[*nbite] = check;
  }
  free(ipiv);
  free(tmp);
}

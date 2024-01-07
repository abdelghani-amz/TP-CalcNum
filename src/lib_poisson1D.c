/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv){
  int ii, jj, kk;
  for (jj=0;jj<(*la);jj++){
    kk = jj*(*lab);
    if (*kv>=0){
      for (ii=0;ii< *kv;ii++){
	AB[kk+ii]=0.0;
      }
    }
    AB[kk+ *kv]=-1.0;
    AB[kk+ *kv+1]=2.0;
    AB[kk+ *kv+2]=-1.0;
  }
  AB[0]=0.0;
  if (*kv == 1) {AB[1]=0;}
  
  AB[(*lab)*(*la)-1]=0.0;
}

void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv){
  int ii, jj, kk;
  for (jj=0;jj<(*la);jj++){
    kk = jj*(*lab);
    if (*kv>=0){
      for (ii=0;ii< *kv;ii++){
	AB[kk+ii]=0.0;
      }
    }
    AB[kk+ *kv]=0.0;
    AB[kk+ *kv+1]=1.0;
    AB[kk+ *kv+2]=0.0;
  }
  AB[1]=0.0;
  AB[(*lab)*(*la)-1]=0.0;
}

void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1){
  int jj;
  RHS[0]= *BC0;
  RHS[(*la)-1]= *BC1;
  for (jj=1;jj<(*la)-1;jj++){
    RHS[jj]=0.0;
  }
}  

void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1){
  int jj;
  double h, DELTA_T;
  DELTA_T=(*BC1)-(*BC0);
  for (jj=0;jj<(*la);jj++){
    EX_SOL[jj] = (*BC0) + X[jj]*DELTA_T;
  }
}  

void set_grid_points_1D(double* x, int* la){
  int jj;
  double h;
  h=1.0/(1.0*((*la)+1));
  for (jj=0;jj<(*la);jj++){
    x[jj]=(jj+1)*h;
  }
}

double relative_forward_error(double* x, double* y, int* la){

  double temp = cblas_ddot(*la, y, 1, y,1);
  temp = sqrt(temp);

  cblas_daxpy(*la, -1.0, x, 1, y, 1);
  double relres = cblas_ddot(*la, y, 1, y,1);
  relres = sqrt(relres);

  return relres / temp;

}

// Factorisation LU
void lu_tridiag(size_t ld, size_t lab, double *A, double *x, double *b)
{

    // Initialization
    for (size_t i = 1; i < ld; i++) {
        A[lab * i + 1] -= A[lab * i] * A[lab * (i - 1) + 2] / A[lab * (i - 1) + 1];
        A[lab * i] /= A[lab * (i - 1) + 1];
    }

    // downward
    x[0] = b[0];
    for (size_t i = 1; i < ld; i++)
        x[i] = b[i] - A[lab * i] * x[i - 1];

    // upward
    b[ld - 1] = x[ld - 1] / A[lab * (ld - 1) + 1];
    for (ssize_t i = ld - 2; i > -1; i--)
        b[i] = (x[i] - A[lab * i + 2] * b[i + 1]) / A[lab * i + 1];

}
int indexABCol(int i, int j, int *lab){
  return j*(*lab)+i;
}

int dgbtrftridiag(int *la, int*n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info){
  return *info;
}

/*You can include any C libraries that you normally use*/
#include "math.h"
#include "mex.h"
// #include "stdlib.h"
// #include "stdio.h"
#define PI 3.141592653589793

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    const mxArray *Ar_p, *Ai_p;
    const mxArray *xr_p, *xi_p;
    double *Ar, *Ai, *xr, *xi, *yr, *yi;
    int i, j, N, MN, M, ind;
    Ar_p = prhs[0];
    Ai_p = prhs[1];
    xr_p = prhs[2];
    xi_p = prhs[3];
    
    Ar = mxGetPr(Ar_p);
    Ai = mxGetPr(Ai_p);
    xr = mxGetPr(xr_p);
    xi = mxGetPr(xi_p);
    
    M = mxGetM(Ar_p);
    N = mxGetN(Ar_p);
    
    plhs[0] = mxCreateDoubleMatrix(M, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(M, 1, mxREAL);
    yr = mxGetPr(plhs[0]);
    yi = mxGetPr(plhs[1]);
    
    for(i = 0; i < M; i++)
    {
        yr[i] = 0;
        yi[i] = 0;
    }
    ind = 0;
    for(j = 0; j < N; j++)
    {
        for(i = 0; i < M; i++)
        {
            
            yr[i] += Ar[ind]*xr[j]-Ai[ind]*xi[j];
            yi[i] += Ar[ind]*xi[j]+Ai[ind]*xr[j];
            ind++;
        }       
    }
}



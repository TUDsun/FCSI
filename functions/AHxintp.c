/*You can include any C libraries that you normally use*/
#include "math.h"
#include "mex.h"   /*--This C library is required*/

#define PI 3.141592653589793

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    mxArray *ftaur_p, *ftaui_p, *EH1_p,*EH2_p, *nL_p, *nL1_p, *dt_p, *prts_p;
    double  *ftaur,   *ftaui,   *EH1, *EH2, *nL, *nL1, *dt, *prts;
    double *outArrayR, *outArrayI, *Wr, *Wi, *outR, *outI;
    double tmp0,tmp1,tmp2,angle;
    int i, j, k, Nr, Nsp, Na, na, ind, Eind, v, r, p, q, K;
    int  Numr,Ni,MnL1,N,ME1,M1,ME2,M2,twoNsp;;
    ftaur_p = prhs[0];
    ftaui_p = prhs[1];
    EH1_p = prhs[2];
    EH2_p = prhs[3];
    nL_p = prhs[4];
    nL1_p = prhs[5];
    dt_p = prhs[6];
    prts_p = prhs[7];
    
    ftaur = mxGetPr(ftaur_p);
    ftaui = mxGetPr(ftaui_p);
    EH1 = mxGetPr(EH1_p);
    EH2 = mxGetPr(EH2_p);
    nL = mxGetPr(nL_p);
    nL1 = mxGetPr(nL1_p);
    dt = mxGetPr(dt_p);
    prts = mxGetPr(prts_p);
    Nr = prts[0];
    Nsp = prts[1];
    Na = prts[2];

    
    Numr =  mxGetM(ftaur_p);
    Ni = Numr/Na;  //Num_fre
    MnL1 = mxGetM(nL1_p);
    N = MnL1/Na;
    ME1 = mxGetM(EH1_p);
    M1 = ME1/Na;
    ME2 =  mxGetM(EH2_p);   //
    M2 = ME2/Na;  //N
    twoNsp = 2*Nsp;
    
    plhs[0] = mxCreateDoubleMatrix(N, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(N, 1, mxREAL);
    outR = mxGetPr(plhs[0]);
    outI = mxGetPr(plhs[1]);
    for(i=0;i<N;i++) {
        outR[i]=0;
        outI[i]=0;
    }
    outArrayR =(double*)malloc((Nr)*sizeof(double));
    outArrayI =(double*)malloc((Nr)*sizeof(double));
    Wr =(double*)malloc((Nr/2)*sizeof(double));
    Wi =(double*)malloc((Nr/2)*sizeof(double));
    for(i=0;i<Nr/2;i++) {
        angle = i*PI*2/Nr;
        Wr[i] = cos(angle);
        Wi[i] = sin(angle);
    }
    i = 0;
    j = 1;
    while(j<Nr) {
        j=j*2;
        i++;
    }
    v = i;
    
    for (na=0;na<Na;na++){
//------------------------------Intp--------------------------------------
        for(i=0;i<Nr;i++) {
            outArrayR[i]=0;
            outArrayI[i]=0;
        }
        
        for (i=0;i<Ni;i++) {
            for (j=1-Nsp;j<=Nsp;j++) {
                ind = nL[i+na*Ni]+j;
                while (ind >= Nr) {
                    ind -= Nr;
                }
                while (ind < 0) {
                    ind += Nr;
                }
                Eind = na*M1+i*twoNsp+j+Nsp-1;
                k = i+na*Ni;
                outArrayR[ind] += ftaur[k]*EH1[Eind];
                outArrayI[ind] += ftaui[k]*EH1[Eind];
            }
        }
///////////////////////////////////倒位序//////////////////////////////////////
        j = 0;
        for( k = 1 ; k<Nr; k++){
            K = Nr/2;
            while(j >= K){//%对j作反方向二进制加法，即j的最高位加1，并向右进位
                j = j-K;
                K = K/2;
            }
            j = j+K;
            if(k < j) {
                tmp0 = outArrayR[k];
                outArrayR[k] = outArrayR[j];
                outArrayR[j] = tmp0;
                tmp0 = outArrayI[k];
                outArrayI[k] = outArrayI[j];
                outArrayI[j] = tmp0;
            }
        }
//////////////////////////////蝶形运算///////////////////////////////////////
        for (i=1;i<=v;i++){
            for (k=1;k<=(1<<(v-i));k++){
                for (j=1;j<=(1<<(i-1));j++){
                    p=(k-1)*(1<<i)+j-1;
                    q=p+(1<<(i-1));
                    r=(j-1)*(1<<(v-i));
                    
                    tmp0=outArrayR[p];
                    tmp1=Wr[r]*outArrayR[q]-Wi[r]*outArrayI[q];
                    tmp2 = Wi[r]*outArrayR[q]+Wr[r]*outArrayI[q];
                    outArrayR[p] = tmp0+tmp1;
                    outArrayR[q] = tmp0-tmp1;
                    
                    tmp0 = outArrayI[p];
                    outArrayI[p] = tmp0+tmp2;
                    outArrayI[q] = tmp0-tmp2;
                }
            }
        }
//////////////////////////////////////////////////////////////////////////
        for(i=0;i<M2;i++){
            ind = na*M2+i;
            outArrayR[i] = outArrayR[i]*EH2[ind]/Nr;
            outArrayI[i] = outArrayI[i]*EH2[ind]/Nr;
        }
        for(i=0;i<N;i++){
            k = i+na*N;
            ind = nL1[k];
            outR[i] += outArrayR[ind]+(outArrayR[ind+1]-outArrayR[ind])*dt[k];
            outI[i] += outArrayI[ind]+(outArrayI[ind+1]-outArrayI[ind])*dt[k];
        }
    }
//-----------------------------------------------------------------------
    free(Wr);
    free(Wi);
    free(outArrayR);
    free(outArrayI);
    
    return;
}
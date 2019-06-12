/*You can include any C libraries that you normally use*/
#include "math.h"
#include "mex.h"   /*--This C library is required*/

#define PI 3.141592653589793

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    const mxArray *ftaur_p, *ftaui_p, *E1_p,*E2_p, *mL_p, *fk_p, *prts_p;
    double  *ftaur, *ftaui, *E1, *E2, *mL, *fk, *prts;
    double *outArrayR, *outArrayI, *Wr, *Wi, *outR, *outI;
    double tmp0,tmp1,tmp2,angle;
    int i, j, k,Na, na, Mr, Msp, ind, Eind,v, r, p, q, K;
    int N,Mfk,Mk,ME1,ME2,M1,M2,twoMsp;
    ftaur_p = prhs[0];
    ftaui_p = prhs[1];
    E1_p = prhs[2];
    E2_p = prhs[3];
    mL_p = prhs[4];
    fk_p = prhs[5];
    prts_p = prhs[6];
    
    ftaur = mxGetPr(ftaur_p);
    ftaui = mxGetPr(ftaui_p);
    E1 = mxGetPr(E1_p);
    E2 = mxGetPr(E2_p);
    mL = mxGetPr(mL_p);
    fk = mxGetPr(fk_p);
    prts = mxGetPr(prts_p);
    Mr = prts[0];
    Msp = prts[1];
    Na = prts[2];

    N =  mxGetM(ftaur_p);
    ME1 = mxGetM(E1_p);
    ME2 =  mxGetM(E2_p);  
    M1 = ME1/Na;
    M2 = ME2/Na;//(kmax+1)
    Mfk = mxGetM(fk_p);
    Mk = Mfk/Na;//Num_fre
    twoMsp = 2*Msp;
    
    
    plhs[0] = mxCreateDoubleMatrix(Mfk, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(Mfk, 1, mxREAL);
    outR = mxGetPr(plhs[0]);
    outI = mxGetPr(plhs[1]);
    outArrayR =(double*)malloc((Mr)*sizeof(double));
    outArrayI =(double*)malloc((Mr)*sizeof(double));
    Wr =(double*)malloc((Mr/2)*sizeof(double));
    Wi =(double*)malloc((Mr/2)*sizeof(double));
    for(i=0;i<Mr/2;i++) {
        angle = -i*PI*2/Mr;
        Wr[i] = cos(angle);
        Wi[i] = sin(angle);
    }
    i = 0;
    j = 1;
    while(j<Mr) {
        j=j*2;
        i++;
    }
    v = i;
    for (na=0;na<Na;na++){
//------------------------------Intp--------------------------------------        
        for(i=0;i<Mr;i++) {
            outArrayR[i]=0;
            outArrayI[i]=0;
        }
        for (i=0;i<N;i++) {
            for (j=1-Msp;j<=Msp;j++) {
                ind = mL[na*N+i]+j;
                while (ind >= Mr) {
                    ind -= Mr;
                }
                while (ind < 0) {
                    ind += Mr;
                }
                Eind = na*M1+i*2*Msp+j+Msp-1;
                outArrayR[ind] += ftaur[i]*E1[Eind];
                outArrayI[ind] += ftaui[i]*E1[Eind];
            }
        }
        
///////////////////////////////////倒位序//////////////////////////////////////
        
        j = 0;
        for( k = 1 ; k<Mr; k++){
            K = Mr/2;
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
///////////////////////////////////////////////////////////////////////////        
        for(i=0;i<Mk;i++){
            ind = na*Mk+i;
            j = fk[ind];
            Eind = j+na*M2;
            outR[ind] = outArrayR[j]*E2[Eind]/Mr;
            outI[ind] = outArrayI[j]*E2[Eind]/Mr;
        }
//-----------------------------------------------------------------------      
    }

    free(Wr);
    free(Wi);
    free(outArrayR);
    free(outArrayI);
    
    return;
}

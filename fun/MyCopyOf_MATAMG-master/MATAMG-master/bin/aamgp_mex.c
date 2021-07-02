/*=========================================================
 * aamgp_mex.c
 *
 * This is a MEX-file for MATAMG.
 *=======================================================*/

#include <math.h>
#include "mex.h"

double getA(double *ia, double *ja, double *aa, int row, int col) {
    double a_i_j = 0.0;
    
    int j;
    for ( j = (int)ia[ row - 1 ]; j < (int)ia[ row ]; j++ ) {
        if ( (int)ja[ j - 1 ] == col ) {
            a_i_j = aa[ j - 1 ];
            break;
        }
    }        
    return a_i_j;
    
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    double *C, *ia, *ja, *aa, *ip,*jp,*ap,*carray,*x;
    int pindex = 1;
    int i,j,jj,k,l;
    int *Ci, *Fi, Ci_index, Fi_index;
    double a_k_j, a_i_k, *temp_Ci, bottom, bottom2;
       
    /* Number of row of matrix A */
    int m;
    
    if(nrhs!=9) {
        mexErrMsgIdAndTxt("MATAMG:aamgp_mex:nrhs", "AAMGP_MEX(C,IA,JA,AA,IP,JP,AP,CARRAY,X");
    }
    
    /* get pointer  */
    C = mxGetPr(prhs[0]);    
    
    ia = mxGetPr(prhs[1]);
    ja = mxGetPr(prhs[2]);
    aa = mxGetPr(prhs[3]);
    
    ip = mxGetPr(prhs[4]);
    jp = mxGetPr(prhs[5]);
    ap = mxGetPr(prhs[6]);
    
    carray = mxGetPr(prhs[7]);
    x = mxGetPr(prhs[8]);
    
    
    m = mxGetN(prhs[0]);

    
    for ( i = 1; i <= m; i++ ){
        if ( C[ i - 1 ] == 1 ) {
            /*mexPrintf("[%d] is C point\n",i);*/
            ip[ pindex - 1 ] = i;
            jp[ pindex - 1 ] = carray[ i - 1 ];
            ap[ pindex - 1 ] = 1.0;
            pindex += 1;            
        }
        else {
            /*mexPrintf("[%d] is F point\n",i);*/
            Ci = mxCalloc((int)ia[ i ] - (int)ia[ i - 1 ],sizeof(int));
            temp_Ci = mxCalloc((int)ia[ i ] - (int)ia[ i - 1 ],sizeof(double));
            Fi = mxCalloc((int)ia[ i ] - (int)ia[ i - 1 ],sizeof(int));
                        
            Ci_index = 1;
            Fi_index = 1;
                      
            bottom2 = 0;
            /* Discrimate Ci,Dsi and Dwi among neighborhood points */
            for ( j = (int)ia[ i - 1 ]; j < (int)ia[ i ]; j++) {
                /*mexPrintf("nbd of %d  = %d\n",i,(int)ja[ j - 1 ]);*/
                if ( C[ (int)ja[ j - 1 ] - 1 ] == 1.0 ) {
                    Ci[ Ci_index - 1 ] = (int)ja[ j - 1 ];
                    ip[ pindex + Ci_index - 2 ] = i;
                    jp[ pindex + Ci_index - 2 ] = carray[ (int)ja[ j - 1 ] - 1 ];
                    ap[ pindex + Ci_index - 2 ] = aa[ j - 1 ];
                    Ci_index += 1;
                    /*mexPrintf("\t[%d] => C_%d\n",(int)ja[ j - 1 ],i);*/
                }
                else if ( i != ja[ j - 1 ] ){
                    Fi[ Fi_index - 1 ] = (int)ja[ j - 1 ];
                    Fi_index += 1;                    
                }
                else {
                    bottom2 += aa[ j - 1 ];
                }
            } /* end loop j */
            
            /* calculate weight P(i,j)*/
            for ( k = 1; k < Fi_index; k++ ) {
                bottom = 0.0;
                a_i_k = getA(ia,ja,aa,i,Fi[ k - 1 ]);
                
                for ( j = 1; j < Ci_index; j++ ) {
                    a_k_j = getA(ia,ja,aa,Fi[ k - 1 ],Ci[ j - 1 ]);
                    if ( a_k_j != 0.0 ) {
                        temp_Ci[ j - 1 ] += a_i_k*a_k_j*x[ Fi[ k - 1 ] - 1 ];
                        bottom += a_k_j * x[ Ci[ j - 1 ] - 1 ];
                    }
                }
                
                for ( j = 1; j < Ci_index; j++ ) {
                    ap[ pindex + j - 2 ] += temp_Ci[j - 1]/bottom;
                    temp_Ci[ j - 1 ] = 0;
                    
                }
            }
            
            for ( j = 0; j < Ci_index; j++) {
                ap[ pindex + j - 1] = - ap[ pindex + j -1 ]/bottom2;                
            }
            pindex += Ci_index - 1;
            mxFree(Ci);
            mxFree(Fi);            
            mxFree(temp_Ci);
        }
    }
    plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
    *mxGetPr(plhs[0])=(double)pindex;
    return;
}

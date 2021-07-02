/*=========================================================
 * aamgp_mex.c
 *
 * This is a MEX-file for MATAMG.
 *=======================================================*/

#include <math.h>
#include "mex.h"
#include "matamg.h"

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
    
    double *C, *ia, *ja, *aa, *ip,*jp,*apr,*carray;
    
    
    int pindex = 1;
    
    /* Number of row of matrix A */
    int m;
    
    int i,j,jj,k,l,isfirst;
    int *Ci, *Fi,  Ci_index, Fi_index;
    double d, *a_i_k;
       
    int percent = 10;   
    
    if(nrhs!=9) {
        mexErrMsgIdAndTxt("MATAMG:ibamgp_mex:nrhs", "IBAMGP_MEX(C,ia,ja,aa,ip,jp,ap,carray,x)");
    }
        
    
    mexPrintf("\t[");
    
    /* get pointer  */
    C = mxGetPr(prhs[0]);    
    
    ia = mxGetPr(prhs[1]);
    ja = mxGetPr(prhs[2]);
    aa = mxGetPr(prhs[3]);
    
    ip = mxGetPr(prhs[4]);
    jp = mxGetPr(prhs[5]);
    apr = mxGetPr(prhs[6]);
    
    carray = mxGetPr(prhs[7]);
    
    
    m = mxGetN(prhs[0]);

    for ( i = 1; i <= m; i++ ){
        if ( (int)(i * 100/m ) == percent ){
            mexPrintf("*");
            percent +=10;
        }
        if ( C[ i - 1 ] == 1 ) {
            ip[ pindex - 1 ] = i;
            jp[ pindex - 1 ] = carray[ i - 1 ];
            apr[ pindex - 1 ] = 1.0;
            pindex += 1;            
        }
        else {
            mxArray *xci, *xfi, *pinvxci,*w;
            
            Ci = mxCalloc((int)ia[ i ] - (int)ia[ i - 1 ],sizeof(int));
            Fi = mxCalloc((int)ia[ i ] - (int)ia[ i - 1 ],sizeof(int));
            a_i_k = mxCalloc((int)ia[ i ] - (int)ia[ i - 1 ],sizeof(double));
                        
            Ci_index = 1;
            Fi_index = 1;
                      
            /* Discrimate Ci,Dsi and Dwi among neighborhood points */
            for ( j = (int)ia[ i - 1 ]; j < (int)ia[ i ]; j++) {
                if ( C[ (int)ja[ j - 1 ] - 1 ] == 1 ) {
                    Ci[ Ci_index - 1 ] = (int)ja[ j - 1 ];
                    Ci_index += 1;
                }
                else if ( i != (int)ja[ j - 1 ] ){
                    Fi[ Fi_index - 1 ] = (int)ja[ j - 1 ];
                    for ( l = (int)ia[ i - 1 ]; l < (int)ia[ i ]; l++ ){
                        if ( ja[ l - 1 ] == ja[ j - 1 ] ) {
                            a_i_k[ Fi_index - 1 ] = aa[ l - 1 ];
                            break;
                        }
                    } /* end loop l */
                        
                    Fi_index += 1;                    
                }
                else {
                    d = aa[ j - 1 ];
                }
            } /* end loop j */
            
                        
            xci = subcol(prhs[8], Ci, Ci_index-1);
            
            xfi = subibamg(prhs[8],a_i_k, Fi, Fi_index-1); 
            
            pinvxci = pinv(xci);
            
            w = times(pinvxci,xfi);
            
            /* calculate weight P(i,j)*/
            for ( j = 1; j < Ci_index; j++ ) {
                ip[ pindex - 1 ] = i;
                jp[ pindex - 1 ] = (int)carray[ Ci[ j - 1 ] - 1 ];
                apr[ pindex - 1 ] = -(getA(ia,ja,aa,i,Ci[ j - 1 ]) + getR(w,j - 1,0))/d;
                pindex += 1;
            } /* end loop j */
            
            
            mxFree(Ci);
            mxFree(Fi);
            mxFree(a_i_k);
            
            mxDestroyArray(xci);
            mxDestroyArray(xfi);
            mxDestroyArray(pinvxci);
            mxDestroyArray(w);
            
        }
        
    }
    
    mexPrintf("]\n");
    plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
    *mxGetPr(plhs[0])=(double)pindex;
    return;
}

/*=========================================================
 * rbamgp_mex.c
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
    
    int i,j,k,l,isfirst;
    int *Ci, Ci_index;
    
    int nofx;
    
    int percent = 10;
    
    if(nrhs!=10) {
        mexErrMsgIdAndTxt("MATAMG:ramgp_mex:nrhs", "RBAMGP_MEX(C,IA,JA,AA,IP,JP,AP,CARRAY,X,R");
    }
      
    
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
    
    
    nofx = mxGetM(prhs[8]);

    
    for ( i = 1; i <= m; i++ ){
        if ( C[ i - 1 ] == 1 ) {
            ip[ pindex - 1 ] = i;
            jp[ pindex - 1 ] = carray[ i - 1 ];
            apr[ pindex - 1 ] = 1.0;            
            pindex += 1;            
        }
        else {
            mxArray *xci,*xi,*w,*pinvxci;
            
            Ci = mxCalloc((int)ia[ i ] - (int)ia[ i - 1 ],sizeof(int));
                        
            Ci_index = 1;
            
            /* Find Ci points */
            for ( j = (int)ia[ i - 1 ]; j < (int)ia[ i ]; j++) {
                if ( C[ (int)ja[ j - 1 ] - 1 ] == 1.0 ) {
                    Ci[ Ci_index - 1 ] = (int)ja[ j - 1 ];
                    Ci_index += 1;
                }
            } /* end loop j */
            
            if ( Ci_index == 1 ) {
                mexPrintf("[%d] has no neighboring C point\n",i);
                for ( j = (int)ia[ i - 1 ]; j < (int)ia[ i ]; j++ ) {
                    for ( k = (int)ia[ (int)ja[ j - 1 ] - 1 ]; k < (int)ia[ (int)ja[ j - 1 ] ]; k++ ){
                        if (C[ (int)ja[ k - 1 ] - 1 ] ) {
                            isfirst = 1;
                            for ( l = 1; l <= Ci_index; l++ ){
                                if ( Ci[ l - 1 ] == (int)ja[ k - 1 ] ) {
                                    isfirst = 0;
                                }
                            }
                            if ( isfirst == 1 ) {
                                Ci[ Ci_index - 1 ] = (int)ja[ k - 1 ];
                                Ci_index += 1;
                            }
                        }
                    }
                }
            }

            xci = subcol(prhs[8],Ci,Ci_index-1);
 
            xi = subrbamg(prhs[8],prhs[9],getA(ia,ja,aa,i,i),0,nofx-1,i-1,i-1);
            
            pinvxci = pinv(xci);
            
            w = times(pinvxci,xi);
            
            
            /* calculate weight P(i,j)*/
            for ( j = 1; j < Ci_index; j++ ) {
                ip[ pindex - 1 ] = i;
                jp[ pindex - 1 ] = (int)carray[ Ci[ j - 1 ] - 1 ];
                
                apr[ pindex - 1 ] = getR(w,j - 1,0);
                pindex += 1;
            } /* end loop j */
            
            mxFree(Ci);
                        
            mxDestroyArray(xci);
            mxDestroyArray(xi);
            mxDestroyArray(pinvxci);
            mxDestroyArray(w);
        }
    }
    
    plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
    *mxGetPr(plhs[0])=(double)pindex;
    return;
}

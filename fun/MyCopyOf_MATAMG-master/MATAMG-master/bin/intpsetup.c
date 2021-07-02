/*=========================================================
 * interp_setup.c
 *
 * This is a MEX-file for MATAMG.
 *=======================================================*/

#include <math.h>
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
        
    double *C, *ia, *carray;
    int m;
    int i;
    int maxnnz = 0;
    
    if(nrhs!=3) {
        mexErrMsgIdAndTxt("MATAMG:interp_setup:nrhs", "INTERP_SETUP(C,IA,CARRAY");
    }
    
    /* get pointer  */
    C = mxGetPr(prhs[0]);    
    ia = mxGetPr(prhs[1]);
    carray = mxGetPr(prhs[2]);
    
    /* Number of row of matrix A */
    m = mxGetN(prhs[0]);

    
    
    /* Measure number of nonzeros of P and define map from coarse to fine grid */
    for ( i = 1; i <= m; i++ ){
        if ( C[ i - 1 ] ) {
            maxnnz += 1;
        }
        else {
            maxnnz += (int)ia[ i ] - (int)ia[ i - 1 ];
        }
        
        /* Find coarse level index for C points */
        if( i == 1 ){
            carray[ i - 1 ] = C[ i - 1 ];
        }
        else {
            carray[ i - 1 ] = carray[ i - 2 ] + C[ i - 1 ];
        }        
    }
    
    plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
    *mxGetPr(plhs[0])=(double)maxnnz;
    return;
}

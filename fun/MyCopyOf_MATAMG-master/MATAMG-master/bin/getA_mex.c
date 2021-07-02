/*=========================================================
 * getA_mex.c
 * Return matrix element
 *
 *
 * get A(row,col)
 *
 * This is a MEX-file for MATAMG.
 *=======================================================*/

#include <math.h>
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    /* check for proper number of arguments*/
    if(nrhs!=5) {
        mexErrMsgIdAndTxt("MATAMG:getA_mex:nrhs", "Usage : a_i_j = getA_mex(IA,JA,AA,ROW,COL)");
    }
    
    double *ia, *ja, *aar, *aai;  /* Compressed row stroage for matirx A */
    int row, col;
    
    /* get pointers to the row pointer array */
    ia = mxGetPr(prhs[0]);
    
    /* get pointers to the column index array */
    ja = mxGetPr(prhs[1]);
    
    /* get pointers to value array */
    aar = mxGetPr(prhs[2]);
    aai = mxGetPi(prhs[2]);
    
    /* get row and column index */
    row  = mxGetScalar(prhs[3]);
    col  = mxGetScalar(prhs[4]);
    
    int j;   
        
    if ( aai != NULL ) {
        plhs[0] = mxCreateDoubleMatrix(1 , 1, mxCOMPLEX);
        *mxGetPr(plhs[0]) = 0.0;
        *mxGetPi(plhs[0]) = 0.0;    
    }
    
    else {
        plhs[0] = mxCreateDoubleMatrix(1 , 1, mxREAL);
        *mxGetPr(plhs[0]) = 0.0;    
    }
    
    for ( j = ia[ row - 1 ]; j < ia[ row ] ; j++ ) {
        if ( (int)ja[ j - 1 ] == col ){
            if ( aai != NULL ){
                *mxGetPr(plhs[0]) = aar[ j - 1 ];
                *mxGetPi(plhs[0]) = aai[ j - 1 ];
            }
            else {
                *mxGetPr(plhs[0]) = aar[ j - 1 ];
            }
            
            break;
        }
    }
    
    return;
}
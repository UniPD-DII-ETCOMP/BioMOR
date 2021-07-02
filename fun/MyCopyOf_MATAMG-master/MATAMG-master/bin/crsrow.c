/*=========================================================
 * crsrow.c
 * Find row index of compressed row stroage
 *
 * The calling syntax is:
 * 
 *
 * This is a MEX-file for MATAMG.
 *=======================================================*/

#include <math.h>
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    
    double *I;
    double *diff;
    double *IA;
    int lengthI, k, m;
    (void) plhs;    /* unused parameter */
    
    /* Print warning if number of right hand side is not three */
    if(nrhs!=3) {
        mexErrMsgIdAndTxt("MATAMG:crsrow:nrhs", "CRSROW(I,DIFF,IA)");
    }
    
    /* get pointer to the row index */
    I = mxGetPr(prhs[0]);
    
    /* get nnz array */
    diff = mxGetPr(prhs[1]);
    
    /* get pointer to ia */
    IA = mxGetPr(prhs[2]);
    
    
    lengthI = mxGetM(prhs[0]);
    m = mxGetM(prhs[1]);
    
    for ( k = 1; k < lengthI + 1; k++ ) {
        diff[ (int)I[ k - 1 ] -1 ] += 1;
    }
    
    IA[0] = 1;
    
    for ( k = 1; k < m + 1; k++ ) {
        IA[ k ]  = IA[ k - 1 ] + diff[ k - 1];  
    }
    return;
}
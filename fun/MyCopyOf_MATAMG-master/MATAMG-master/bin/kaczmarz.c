/*=========================================================
 * kaczmarz.c
 * Karczmarx method for symmetric matrix
 *
 *
 * iterates Ax=b
 *
 * This is a MEX-file for MATAMG.
 *=======================================================*/

#include <math.h>
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    int maxiter;            /* Maximum number of iteration */
    int m;
    double *ia, *ja, *aar, *aai;    /* Compressed row stroage for matirx A */
    double *vr, *vi;                /* initial vector */
    double *fr, *fi;                /* RHS vector */
    double w;                       /* Iteration parameter */
    double valr, vali;              /* temporary value for summation */
    double *rownorm_r, *rownorm_i;
    
    int i, j, k;
    double topr, topi;
    
    (void) plhs;
    
    /* check for proper number of arguments*/
    if(nrhs < 6 || nrhs > 7) {
        mexErrMsgIdAndTxt("MATAMG:kaczmarz:nrhs", "Usage : V = KACZMARZ(IA,JA,AA,V,F,W,MAXITER,ROWNORM)");
    }
    
    /* get pointers to the row pointer array */
    ia = mxGetPr(prhs[0]);
    
    /* get pointers to the column index array */
    ja = mxGetPr(prhs[1]);
    
    /* get pointers to value array */
    aar = mxGetPr(prhs[2]);
    aai = mxGetPi(prhs[2]);
    
    /* get the inital vector */
    vr  = mxGetPr(prhs[3]);
    vi  = mxGetPi(prhs[3]);
    
    /* mexPrintf("Pointer : %f\n", vi);*/
    /* get RHS vector */
    fr  = mxGetPr(prhs[4]);
    fi  = mxGetPi(prhs[4]);
    
    
    /* get maximum number of iteration sweeps*/
    maxiter = mxGetScalar(prhs[5]);
    
    /* get the size of matrix */
    m = mxGetM(prhs[0]) -1;
    
    
    
    if ( nrhs == 6 ){
        rownorm_r = mxCalloc(m,sizeof(double));
        if (aai == NULL) {
            for ( i = 1; i <= m; i++ ) {
                rownorm_r[ i - 1 ] = 0.0;
                for ( j = (int)ia[ i - 1 ]; j < (int)ia[ i ]; j++ ) {
                    rownorm_r[ i - 1 ] += aar[ j - 1 ]*aar[ j - 1 ];
                }
            }
        }
        else {
            /* Complex part */
        }
    }
    
    
    
    
    
    if ( aai == NULL )      /* Check whether matrix is complex or not */ {   /* Real matrix and real vector */
        if ( maxiter > 0 ) {
            for ( k = 1; k <= maxiter; k++ ) {
                for ( i = 1; i < m + 1 ; i++ ) {
                    topr = -fr[ i - 1 ]; 
                    
                    for ( j = (int)ia[ i - 1 ]; j < (int)ia[ i ] ; j++ ) {
                        topr += aar[ j - 1 ] * vr[ (int)ja[ j - 1 ] - 1 ];                        
                    } /* end loop j */
                    
                    for ( j = (int)ia[ i - 1 ]; j < (int)ia[ i ] ; j++ ) {
                        vr[ (int)ja[ j - 1 ] - 1 ] -= topr * aar[ j - 1 ] / rownorm_r[ i - 1 ];
                    } /* end loop j */;
                    
                } /* end loop i */
            } /* end loop k */
        }
    }
    else {   /* Complex matrix */
        
    }
    
    if (nrhs == 6 ) {
        if ( aai == NULL ) {
            mxFree(rownorm_r);
        }
        else {
        /* complex */   
        }
    }
            
            
            
    return;
}



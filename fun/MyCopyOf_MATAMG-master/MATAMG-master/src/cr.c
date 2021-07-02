/* gs.c */


#include <math.h>
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    int maxiter;            /* Maximum number of iteration */
    int sizeofA;
    double *ia, *ja, *aar, *aai;  /* Compressed row stroage for matirx A */
    double *vr, *vi;              /* initial vector */
    double *fr, *fi;              /* RHS vector */
    double w;               /* Iteration parameter */
    double *C;              /* Array for coarse grid */
    double valr, vali;
    
    int i, j, k;
    double dr, di;
    
    (void) plhs;
    
    /* check for proper number of arguments*/
    if(nrhs!=8) {
        mexErrMsgIdAndTxt("MATAMG:cr:nrhs", "Usage : V = CR(IA,JA,AA,V,F,W,MAXITER,C)");
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
    
    /* get RHS vector */
    fr  = mxGetPr(prhs[4]);
    fi  = mxGetPi(prhs[4]);
    
    /* get iteration parameter */
    w  = mxGetScalar(prhs[5]);
    
    /* get maximum number of iteration sweeps*/
    maxiter = mxGetScalar(prhs[6]);
    
    /* get array C */
    C = mxGetPr(prhs[7]);
    
    /* get the size of matrix */
    sizeofA = mxGetM(prhs[0]) -1;
    
    /* check the size of matrix A and C */
    if (sizeofA != mxGetN(prhs[7]) ) {
        mexPrintf("C : %d, A : %d \n", mxGetN(prhs[7]), sizeofA);
        mexErrMsgIdAndTxt("MATAMG:cr:C", "Size of C doesn't match that of A");
    }
    
    
    if ( aai == NULL ) {
        for ( i = 1; i < sizeofA + 1; i++ ) {
            if (C[i-1]) {
                vr[ i - 1 ] = 0.0;
            }
        }
    }
    else {
        for ( i = 1; i < sizeofA + 1; i++ ) {
            if (C[i-1]) {
                vr[ i - 1 ] = 0.0;
                vi[ i - 1 ] = 0.0;
            }
        }
    }
    if ( aai == NULL )      /* Real valued matrix */ {
        if ( maxiter > 0 ) {
            for ( k = 1; k <= maxiter; k++ ) {
                /* F - Relaxation */
                for ( i = 1; i < sizeofA + 1 ; i++ ) {
                    if ( !C[ i - 1] ) {
                        valr = 0.0; vali = 0.0;
                        dr = 0.0; di = 0.0;
                        for ( j = ia[ i - 1 ]; j < ia[ i ] ; j++ ) {
                            if( i != (int)ja[ j - 1 ]  ) {
                                valr = valr + aar[ j - 1 ] * vr[ (int)ja[ j - 1 ] -1 ];
                            }
                            else {
                                dr = aar[ j - 1 ];
                            }
                            
                        } /* end loop j */
                        if (dr == 0.0)
                            mexErrMsgIdAndTxt("MATAMG:gs:diag", "There is zero diagonal");
                        vr[ i - 1 ] = w * ( fr[ i - 1 ] - valr )/dr + ( 1 - w ) * vr[ i -1 ];
                    }
                } /* end loop i */
            } /* end loop k */
        }
    }
    
    else {   /* Complex matrix */
        if ( vi == NULL ) {
            mexErrMsgIdAndTxt("MATAMG:gs:complex", "Initial vector must be complex");
        }
        else {
            if ( fi == NULL ) {
                if ( maxiter > 0 ) {
                    for ( k = 1; k <= maxiter; k++ ) {
                        for ( i = 1; i < sizeofA + 1 ; i++ ) {
                            /* F - Relaxation */
                            if ( !C[ i - 1] ) {
                                valr = 0.0;
                                vali = 0.0;
                                dr = 0.0;
                                di = 0.0;
                                for ( j = ia[ i - 1 ]; j < ia[ i ] ; j++ ) {
                                    if( i != (int)ja[ j - 1 ]  ) {
                                        valr = valr + aar[ j - 1 ] * vr[ (int)ja[ j - 1 ] -1 ] - aai[ j - 1 ] * vi[ (int)ja[ j - 1 ] -1 ];
                                        vali = vali + aar[ j - 1 ] * vi[ (int)ja[ j - 1 ] -1 ] + aai[ j - 1 ] * vr[ (int)ja[ j - 1 ] -1 ];
                                    }
                                    else {
                                        dr = aar[ j - 1 ];
                                        di = aai[ j - 1 ];
                                    }
                                } /* end loop j */
                                if (dr == 0.0 & di == 0.0 )
                                    mexErrMsgIdAndTxt("MATAMG:gs:diag", "There is zero diagonal");
                                vr[ i - 1 ] = w * (( fr[ i - 1 ] - valr)*dr -  vali*di)/(dr*dr+di*di) + (1-w)*vr[ i - 1 ];
                                vi[ i - 1 ] = w * (( valr - fr[ i - 1 ] )*di - vali*dr)/(dr*dr+di*di) + (1-w)*vi[ i - 1 ];
                            }
                        } /* end loop i */
                    } /* end loop k */
                }
            }
            else {
                if ( maxiter > 0 ) {
                    for ( k = 1; k <= maxiter; k++ ) {
                        for ( i = 1; i < sizeofA + 1 ; i++ ) {
                            if ( !C[ i - 1 ] ){
                                
                                valr = 0.0;
                                vali = 0.0;
                                dr = 0.0;
                                di = 0.0;
                                for ( j = ia[ i - 1 ]; j < ia[ i ] ; j++ ) {
                                    if( i != (int)ja[ j - 1 ]  ) {
                                        valr = valr + aar[ j - 1 ] * vr[ (int)ja[ j - 1 ] -1 ] - aai[ j - 1 ] * vi[ (int)ja[ j - 1 ] -1 ];
                                        vali = vali + aar[ j - 1 ] * vi[ (int)ja[ j - 1 ] -1 ] + aai[ j - 1 ] * vr[ (int)ja[ j - 1 ] -1 ];
                                    }
                                    else {
                                        dr = aar[ j - 1 ];
                                        di = aai[ j - 1 ];
                                    }
                                } /* end loop j */
                                if (dr == 0.0 & di == 0.0 )
                                    mexErrMsgIdAndTxt("MATAMG:gs:diag", "There is zero diagonal");
                                vr[ i - 1 ] = w * (( fr[ i - 1 ] - valr)*dr + ( fi[ i - 1 ] - vali)*di)/(dr*dr+di*di) + (1-w)*vr[ i - 1 ];
                                vi[ i - 1 ] = w * (( valr - fr[ i - 1 ] )*di + ( fi[ i - 1 ] - vali)*dr)/(dr*dr+di*di) + (1-w)*vi[ i - 1 ];
                            }
                        }/* end loop i*/
                    }
                }
                
            } /* end else */
        }
    }
    return;
}
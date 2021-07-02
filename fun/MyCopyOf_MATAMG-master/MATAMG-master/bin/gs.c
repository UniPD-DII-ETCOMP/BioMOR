/*=========================================================
 * gs.c
 * Gauss-Seidel method
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
    int sizeofA;
    double *ia, *ja, *aar, *aai;    /* Compressed row stroage for matirx A */
    double *vr, *vi;                /* initial vector */
    double *fr, *fi;                /* RHS vector */
    double w;                       /* Iteration parameter */
    double valr, vali;              /* temporary value for summation */
    
    int i,im1, j,jm1, k, jaj, iai, iaim1;
    double dr, di, bottom;
    
    (void) plhs;
    
    /* check for proper number of arguments*/
    if(nrhs!=7) {
        mexErrMsgIdAndTxt("MATAMG:gs:nrhs", "Usage : V = GS(IA,JA,AA,V,F,W,MAXITER)");
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
    
    /* get iteration parameter */
    w  = mxGetScalar(prhs[5]);
    
    /* get maximum number of iteration sweeps*/
    maxiter = mxGetScalar(prhs[6]);
    
    /* get the size of matrix */
    sizeofA = mxGetM(prhs[0]) -1;
    
    
    if ( aai == NULL )      /* Check whether matrix is complex or not */ {   /* Real matrix and real vector */
        if ( maxiter > 0 ) {
            for ( k = 1; k <= maxiter; k++ ) {
                iaim1 = ia[ 0 ];
                for ( i = 1; i < sizeofA + 1 ; i++ ) {
                    iai = ia[ i ];
                    valr = 0.0; 
                    dr = 0.0; 
                    im1 = i - 1;
                    for ( j = iaim1; j < iai; j++ ) {
                        jm1 = j - 1;
                        jaj = (int)ja[ jm1 ];
                        if( i != jaj  ) {
                            valr = valr + aar[ jm1 ] * vr[ jaj -1 ];
                        }
                        else {
                            dr = aar[ jm1 ];
                        }
                        
                    } /* end loop j */
                    if (dr == 0.0)
                        mexErrMsgIdAndTxt("MATAMG:gs:diag", "There is zero diagonal");
                    vr[ im1 ] = w * ( fr[ im1 ] - valr )/dr + ( 1 - w ) * vr[ im1 ];
                    iaim1 = iai;
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
                        iaim1 = ia[ 0 ];
                        for ( i = 1; i < sizeofA + 1 ; i++ ) {
                            iai = ia[ i ];
                            valr = 0.0;
                            vali = 0.0;
                            dr = 0.0;
                            di = 0.0;
                            im1 = i - 1;
                            for ( j = iaim1; j < iai; j++ ) {
                                jaj = (int)ja[ j - 1 ];
                                jm1 = j - 1;
                                if( i != jaj  ) {
                                    valr = valr + aar[ jm1 ] * vr[ jaj -1 ] - aai[ jm1 ] * vi[ jaj -1 ];
                                    vali = vali + aar[ jm1 ] * vi[ jaj -1 ] + aai[ jm1 ] * vr[ jaj -1 ];
                                }
                                else {
                                    dr = aar[ jm1 ];
                                    di = aai[ jm1 ];
                                }
                            } /* end loop j */
                            if (dr == 0.0 & di == 0.0 )
                                mexErrMsgIdAndTxt("MATAMG:gs:diag", "There is zero diagonal");
                            bottom = dr*dr + di*di;
                            vr[ im1 ] = w * (( fr[ im1 ] - valr)*dr -  vali*di)/bottom + (1-w)*vr[ im1 ];
                            vi[ im1 ] = w * (( valr - fr[ im1 ] )*di - vali*dr)/bottom + (1-w)*vi[ im1 ];
                            iaim1 = iai;
                        }                 
                    }
                }
            }
            else {
                if ( maxiter > 0 ) {
                    for ( k = 1; k <= maxiter; k++ ) {
                        iaim1 = ia[ 0 ];
                        for ( i = 1; i < sizeofA + 1 ; i++ ) {
                            iai = ia[ i ];
                            valr = 0.0;
                            vali = 0.0;
                            dr = 0.0;
                            di = 0.0;
                            im1 = i - 1;
                            for ( j = iaim1; j < iai ; j++ ) {
                                jaj = (int)ja[ j - 1 ];
                                jm1 = j - 1;
                                if( i != jaj  ) {
                                    valr = valr + aar[ jm1 ] * vr[ jaj -1 ] - aai[ jm1 ] * vi[ jaj -1 ];
                                    vali = vali + aar[ jm1 ] * vi[ jaj -1 ] + aai[ jm1 ] * vr[ jaj -1 ];
                                }
                                else {
                                    dr = aar[ jm1 ];
                                    di = aai[ jm1 ];
                                }
                            } /* end loop j */
                            if (dr == 0.0 & di == 0.0 )
                                mexErrMsgIdAndTxt("MATAMG:gs:diag", "There is zero diagonal");
                            bottom = dr*dr + di*di;
                            vr[ im1 ] = w * (( fr[ im1 ] - valr)*dr + ( fi[ im1 ] - vali)*di)/bottom + (1-w)*vr[ im1 ];
                            vi[ im1 ] = w * (( valr - fr[ im1 ] )*di + ( fi[ im1 ] - vali)*dr)/bottom + (1-w)*vi[ im1 ];
                            iaim1 = iai;
                        }
                    }
                }
                
            }
        }
    }
    return;
}



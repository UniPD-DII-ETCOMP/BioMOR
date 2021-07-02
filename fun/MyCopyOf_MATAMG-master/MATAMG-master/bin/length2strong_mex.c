/*=========================================================
 * length2strong_mex.c
 *
 * This is a MEX-file for MATAMG.
 *=======================================================*/

#include <math.h>
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double *ia, *ja, *C;
    double *is2, *js2, *as2;
    double *inbd, *jnbd, *anbd;
    int m;
    int path;
    int *npath;
    int nbdindex, nbdsize;
    int *len2nbd, len2start;
    int index = 0;
    int i, j, k, jaj, jak,l;
    int isinnbd;
    
    if(nrhs!=4) {
        mexErrMsgIdAndTxt("MATAMG:length2strong:nrhs", "LENGTH2STRONG_MEX(IA,JA,IS,JS,AS,C,PATH)");
    }
    
    ia = mxGetPr(prhs[0]);
    ja = mxGetPr(prhs[1]);
        
    C = mxGetPr(prhs[2]);
    path =(int)mxGetScalar(prhs[3]);
    m = mxGetN(prhs[2]);
    npath = (int*)calloc(m,sizeof(int));
    int maxnnzrow = 0;
    
    for ( i = 1; i <= m; i++ ) {
        if ((int)(ia[ i ] - ia[ i - 1 ]) > maxnnzrow )
            maxnnzrow = (int)(ia[ i ] - ia[ i - 1 ]);
    }
    
    plhs[0] = mxCreateDoubleMatrix(maxnnzrow*maxnnzrow*m,1,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(maxnnzrow*maxnnzrow*m,1,mxREAL);
    plhs[2] = mxCreateDoubleMatrix(maxnnzrow*maxnnzrow*m,1,mxREAL);
    plhs[3] = mxCreateDoubleMatrix(1,1,mxREAL);
    plhs[4] = mxCreateDoubleMatrix(m+1,1,mxREAL);
    plhs[5] = mxCreateDoubleMatrix(maxnnzrow*maxnnzrow*m,1,mxREAL);
    len2nbd = (int*)calloc(maxnnzrow*maxnnzrow,sizeof(int));
    
    is2 = mxGetPr(plhs[0]);
    js2 = mxGetPr(plhs[1]);
    as2 = mxGetPr(plhs[2]);
    inbd = mxGetPr(plhs[4]);
    jnbd = mxGetPr(plhs[5]);
    
    nbdindex = 0;    
    /* Compute length 2 connections */
    inbd[ 0 ] = 1;
    for ( i = 1; i <= m; i++ ) {
        nbdsize = 0;
        /* Check length 1 connections */
        for (j = (int)ia[ i - 1 ]; j < (int)ia[ i ]; j++){
            jaj = (int)ja[ j - 1];
            if (jaj!=i){
                len2nbd[nbdsize] = (int)ja[ j - 1 ];
                nbdsize++;
            }
        }
        len2start = nbdsize;
        /* Check length 2 connections */
        for (j = (int)ia[ i - 1 ]; j < (int)ia[ i ]; j++){
            jaj = (int)ja[ j - 1 ];
            if ( jaj != i ){
                for ( k = (int)ia[ jaj - 1 ]; k < (int)ia[ jaj ] ; k++ ){
                    isinnbd = 0;
                    jak =  (int)ja[ k - 1 ];
                    if (jak != i && jak != jaj){
                        for ( l = 0; l < nbdsize; l++ ) {
                            if ( len2nbd[ l ] == jak ){
                                isinnbd = 1;
                                break;
                            }
                        }
                        if (!isinnbd){
                            jnbd[ nbdindex ] = jak;
                            len2nbd[ nbdsize ] = jak;
                            nbdsize++;
                            nbdindex++;
                        }
                    }
                }
            }
        }
        inbd[ i  ] = inbd[ i - 1 ] + nbdsize - len2start;
    }
    
    index    = 0;
    for ( i = 1; i <= m; i++ ) {
        if ( C[ i - 1 ] ){
             for ( j = (int)ia[ i - 1 ]; j < (int)ia[ i ]; j++ ) {
                jaj = (int)ja[ j - 1 ];
                if ( !C[ jaj - 1 ] ) {
                    /*mexPrintf("  %d is a F point nbd of %d\n", jaj, i);*/
                    for ( k = (int)ia [ jaj - 1 ]; k < (int)ia[ jaj ] ; k++ ){
                        jak = (int)ja[ k - 1 ];
                        isinnbd = 0;
                        for ( l = inbd[ i - 1 ]; l < inbd[ i ];  l++ ){
                            if (jnbd[ l - 1 ] == jak){
                                isinnbd = 1;
                                break;
                            }
                        }
                        
                        if (isinnbd) {
                            if ( C [ jak - 1 ] && jak != i){
                                npath[ jak - 1 ]++;
                            }
                        }
                    }                        
                }
            }
            
            for ( j = (int)ia[ i - 1 ]; j < (int)ia[ i ]; j++ ) {
                jaj = (int)ja[ j - 1 ];
                if ( !C[ jaj - 1 ] ) {
                    for ( k = (int)ia [ jaj - 1 ]; k < (int)ia[ jaj ] ; k++ ){
                        jak = (int)ja[ k - 1 ];
                        if ( C [ jak - 1 ] && jak != i && npath[ jak - 1 ] > 0  ){
                            if ( npath[ jak - 1 ] >= path){
                                is2[ index ] = i;
                                js2[ index ] = jak;
                                as2[ index ] = 1;
                                index++;
                            }
                            npath[ jak - 1 ] = 0.0;
                            
                        }
                    }                        
                }
            }
        }        
    }    
    free(npath);
    free(len2nbd);
    *mxGetPr(plhs[3])=(double)index;
    return;
}
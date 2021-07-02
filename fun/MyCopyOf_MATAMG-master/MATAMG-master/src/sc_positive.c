/* sc_positive.c */

#include <stdlib.h>
#include <math.h>
#include "mex.h"

#define ABS( x ) ( x >= 0.0 ? x : -x )
#define CABS( x, y ) sqrt(x*x+y*y)


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    
    double *ia, *ja, *aar, *aai;   /* array for CRS form */
    double *is;
    int *js;      /* CRS for strong connection */
    double *it, *jt;
    double *newis, *newjs, *newas;
    double *ns, *nst;
    double theta;           /* threshold */
    int sizeofA;
    int index = 0;
    double maxval;
    int i, j, im1, jm1, jajm1, isindex, jtindex;
    double absaarj,aarj, aaij;
    int **jjt;
    
    /* check for proper number of arguments*/
    if(nrhs!=4) {
        mexErrMsgIdAndTxt("MATAMG:strong_connection:nrhs", "Usage : S = STRONG_CONNECTOINS(IA,JA,AA,THETA)");
    }
    
    /* allocate js */
    js = mxCalloc(mxGetM(prhs[1]), sizeof(int));
    
    /* get pointers to the row pointer array */
    ia = mxGetPr(prhs[0]);
    
    /* get pointers to the column index array */
    ja = mxGetPr(prhs[1]);
    
    /* get pointers to value array */
    aar = mxGetPr(prhs[2]);
    aai = mxGetPi(prhs[2]);
    
    /* get coarsening thereshold */
    theta  = mxGetScalar(prhs[3]);
    
    /* get the size of matrix */
    sizeofA = mxGetM(prhs[0]) -1;
    
    plhs[0] = mxCreateDoubleMatrix(sizeofA , 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(sizeofA , 1, mxREAL);
    
    plhs[2] = mxCreateDoubleMatrix(sizeofA + 1 , 1, mxREAL);
    
    /* number of strong dependencies */
    ns =  mxGetPr(plhs[0]);
    /* number of strong influences */
    nst = mxGetPr(plhs[1]);
    
    is = mxGetPr(plhs[2]);
    
    /* column index array of S transpose */
    jjt = malloc(sizeofA * sizeof(int *));
    
    if (jjt == NULL)
        mexErrMsgTxt("out of memory");
    
    for ( i = 1; i <= sizeofA; i++ ) {
        im1 = i - 1;
        jjt[ im1 ] = malloc((int)(ia[ i ] - ia[ im1 ])*sizeof(int));
        nst[ im1 ] = 0;
        if (jjt[ im1 ] == NULL)
            mexErrMsgTxt("out of memory");
    }
    
    is [ 0 ] = 1;
    
    if ( aai == NULL ) {
        for ( i = 1; i <= sizeofA; i++ ) {
            maxval = 0.0;
            im1 = i - 1;
            
            for ( j = ia[ im1 ]; j < ia[ i ]; j++ ) {
                jm1 = j - 1;
                absaarj = ABS(aar[ jm1 ]);
                if ( i != (int)ja[ jm1 ] && maxval < absaarj ) {
                    maxval = absaarj;
                }
            }
            isindex = 0;
            for ( j = ia[ im1 ]; j < ia[ i ]; j ++ ) {
                jm1 = j - 1;
                absaarj = ABS(aar[ jm1 ]);
                if ( i !=(int)ja[ jm1 ] && absaarj > theta*maxval ) {
                    jajm1 = (int)ja[ jm1 ] - 1;
                    js[ index ] = ja[ jm1 ];
                    ns[ i - 1 ]++;
                    jjt[ jajm1 ][ (int)nst[ jajm1 ] ] = i;
                    nst[ jajm1 ]++;
                    index++;
                    isindex++;
                }
            }
            is[ i ] = is[ i - 1 ] + isindex;
        }
    }
    else {
        for ( i = 1; i <= sizeofA; i++ ) {
            maxval = 0.0;
            im1 = i - 1;
            for ( j = ia[ im1 ]; j < ia[ i ]; j++ ) {
                jm1 = j - 1;
                if ( i != (int)ja[ jm1 ] && maxval < aar[ jm1]*aar[ jm1]+aai[ jm1]*aai[ jm1] ) {
                    maxval = aar[ jm1]*aar[ jm1]+aai[ jm1]*aai[ jm1];
                }
            }
            
            isindex = 0;
            for ( j = ia[ i - 1 ]; j < ia[ i ]; j ++ ) {
                jm1 = j - 1;
                if ( i !=(int)ja[ jm1 ] && aar[ jm1]*aar[ jm1]+aai[ jm1]*aai[ jm1]> theta*maxval ) {
                    is[ index ] = i;
                    js[ index ] = ja[ jm1];
                    ns[ i - 1 ]++;
                    jjt[ jajm1 ][ (int)nst[ jajm1 ] ] = i;
                    nst[ jajm1 ]++;
                    index++;
                    isindex++;
                }
            }
            is[ i ] = is[ i - 1 ] + isindex;
        }
    }
    
    plhs[3] = mxCreateDoubleMatrix(index , 1, mxREAL);
    
    plhs[4] = mxCreateDoubleMatrix(sizeofA + 1 , 1, mxREAL);
    plhs[5] = mxCreateDoubleMatrix(index , 1, mxREAL);
    
    newjs = mxGetPr(plhs[3]);
    
    it = mxGetPr(plhs[4]);
    jt = mxGetPr(plhs[5]);
    
    for( i = 0; i < index; i++ ) {
        newjs[ i ] = js[ i ];
    }
    
    it[ 0 ] = 1;
    jtindex = 0;
    for ( i = 0; i < sizeofA; i++){
        for ( j = 0; j < (int)nst[ i ]; j++){
            jt[ jtindex] = jjt[ i ][ j ];
            jtindex++;
        }
        it[ i + 1 ] = it[ i ] + (int)nst[ i ];
    }
    
    for (i = 0; i < sizeofA; i++){
        free(jjt[ i ]);
    }
    
    free(jjt);
    mxFree(js);
}

/* TODO Dynamic threshold */
/*else if (theta == 0.0) {
 * double top;
 * double bottom;
 * double epsilon;
 * if ( aai == NULL ) {
 * double absaar;
 * for ( i = 1; i <= sizeofA; i++ ) {
 * maxval = 0.0;
 * top = 0.0;
 * bottom = 0.0;
 * for ( j = ia[ i - 1 ]; j < ia[ i ]; j++ ) {
 * if ( i != (int)ja[ j - 1 ] ){
 * absaar = ABS(aar[ j - 1 ]);
 * top += absaar*absaar;
 * bottom += absaar;
 * if ( maxval < absaar ) {
 * maxval = absaar;
 * }
 * }
 * }
 * epsilon = top/(maxval*bottom);
 *
 * for ( j = ia[ i - 1 ]; j < ia[ i ]; j ++ ) {
 * if ( i !=(int)ja[ j - 1 ] && ABS(aar[ j - 1 ]) >= epsilon*maxval ) {
 * is[ index ] = i;
 * js[ index ] = ja[ j - 1];
 * as[ index ] = 1;
 * index++;
 * }
 * }
 * }
 * }
 *
 * else {
 * for ( i = 1; i <= sizeofA; i++ ) {
 * maxval = 0.0;
 * for ( j = ia[ i - 1 ]; j < ia[ i ]; j++ ) {
 * if ( i != (int)ja[ j - 1 ] && maxval < aar[ j - 1]*aar[ j - 1]+aai[ j - 1]*aai[ j - 1] ) {
 * maxval = aar[ j - 1]*aar[ j - 1]+aai[ j - 1]*aai[ j - 1];
 * }
 * }
 * for ( j = ia[ i - 1 ]; j < ia[ i ]; j ++ ) {
 * if ( i !=(int)ja[ j - 1 ] && aar[ j - 1]*aar[ j - 1]+aai[ j - 1]*aai[ j - 1]> theta*maxval ) {
 * is[ index ] = i;
 * js[ index ] = ja[ j - 1];
 * as[ index ] = 1;
 * index++;
 * }
 * }
 * }
 * }
 * }
 */
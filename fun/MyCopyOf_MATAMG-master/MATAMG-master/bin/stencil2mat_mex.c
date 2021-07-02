/*=========================================================
 * stencil2mat_mex.c
 *
 * length(stencil) : 3(1D) [0]C [1]W [2]E
 *                 : 5(2D) [0]C [1]W [2]E [3]S [4]N
 *                 : 9(2D) [0]C [1]W [2]E [3]SW [4]S [5]SE [6]NW [7]N [8]NE
 *                 : 7(3D) [0]C [1]W [2]E [3]S [4]N [5]D [6]U
 * This is a MEX-file for MATAMG.
 *=======================================================*/

#include <math.h>
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    
    double *stencilr, *stencili;
    int m,n,size;
    double *ia, *ja, *aar, *aai;
    int i, j;
    
    if(nrhs!=2) {
        mexErrMsgIdAndTxt("MATAMG:stencil2mat:nrhs", "STENCIL2MAT(STENCIL,N)");
    }
    
    /* get pointer  */
    stencilr = mxGetPr(prhs[0]);
    stencili = mxGetPi(prhs[0]);
    
    
    m = (mxGetM(prhs[0])>mxGetN(prhs[0])?mxGetM(prhs[0]):mxGetN(prhs[0]));
    n = (mxGetM(prhs[0])<mxGetN(prhs[0])?mxGetM(prhs[0]):mxGetN(prhs[0]));
    
    size = (int)(*mxGetPr(prhs[1]));
    if ( n > 1 ) {
        mexErrMsgTxt("Stecil must be a vector");
    }
    
    
    if (stencili == NULL){
        switch ( m )
        {
            case ( 3 ):  {
                plhs[0] = mxCreateDoubleMatrix(3*size-5,1,mxREAL);
                plhs[1] = mxCreateDoubleMatrix(3*size-5,1,mxREAL);
                plhs[2] = mxCreateDoubleMatrix(3*size-5,1,mxREAL);
                
                ia = mxGetPr(plhs[0]);
                ja = mxGetPr(plhs[1]);
                aar = mxGetPr(plhs[2]);
                
                j = 0;
                for ( i = 0; i < size -1; i++ ){
                    ia[j] = i+1;
                    ja[j] = i+1;
                    aar[j] = stencilr[0];
                    j++;
                    
                    if ( i - 1 >= 0 ) {
                        ia[j] = i+1;
                        ja[j] = i;
                        aar[j] = stencilr[1];
                        j++;
                    }
                    
                    if ( i < size - 2 ) {
                        ia[j] = i+1;
                        ja[j] = i+2;
                        aar[j] = stencilr[2];
                        j++;
                    }
                }
                break;
            }
            
            case ( 5 ):  {
                plhs[0] = mxCreateDoubleMatrix(5*size*size-14*size+9,1,mxREAL);
                plhs[1] = mxCreateDoubleMatrix(5*size*size-14*size+9,1,mxREAL);
                plhs[2] = mxCreateDoubleMatrix(5*size*size-14*size+9,1,mxREAL);
                
                ia = mxGetPr(plhs[0]);
                ja = mxGetPr(plhs[1]);
                aar = mxGetPr(plhs[2]);
                
                j = 0;
                for ( i = 0; i < (size -1)*(size -1); i++ ){
                    /* Center value */
                    ia[j] = i+1;
                    ja[j] = i+1;
                    aar[j] = stencilr[0];
                    j++;
                    
                    /* Left value */
                    if ( (i+1)%(size-1)!= 1 ) {
                        ia[j] = i+1;
                        ja[j] = i;
                        aar[j] = stencilr[1];
                        j++;
                    }
                    
                    /* Right value */
                    if ( (i+1)%(size-1)!= 0 ) {
                        ia[j] = i+1;
                        ja[j] = i+2;
                        aar[j] = stencilr[2];
                        j++;
                    }
                    
                    /* Down value */
                    if ( i >= size-1 ) {
                        ia[j] = i+1;
                        ja[j] = i+1-(size-1);
                        aar[j] = stencilr[3];
                        j++;
                    }
                    
                    /* Up value */
                    if ( i < (size-1 )*(size-2) ) {
                        ia[j] = i+1;
                        ja[j] = i+1+ (size-1);
                        aar[j] = stencilr[4];
                        j++;
                    }
                }
                break;
            } /* end case (5) */
            
            case ( 9 ):  {
                plhs[0] = mxCreateDoubleMatrix(3*(3*size-1)*(size-3)+16,1,mxREAL);
                plhs[1] = mxCreateDoubleMatrix(3*(3*size-1)*(size-3)+16,1,mxREAL);
                plhs[2] = mxCreateDoubleMatrix(3*(3*size-1)*(size-3)+16,1,mxREAL);
                
                ia = mxGetPr(plhs[0]);
                ja = mxGetPr(plhs[1]);
                aar = mxGetPr(plhs[2]);
                
                j = 0;
                for ( i = 0; i < (size -1)*(size -1); i++ ){
                    /* Center value */
                    ia[j] = i+1;
                    ja[j] = i+1;
                    aar[j] = stencilr[0];
                    j++;
                    
                    /* Left value */
                    if ( (i+1)%(size-1)!= 1 ) {
                        ia[j] = i+1;
                        ja[j] = i;
                        aar[j] = stencilr[1];
                        j++;
                    }
                    
                    /* Right value */
                    if ( (i+1)%(size-1)!= 0 ) {
                        ia[j] = i+1;
                        ja[j] = i+2;
                        aar[j] = stencilr[2];
                        j++;
                    }
                    
                    /* SW value */
                    if (( i >= size-1 ) && ( (i+1)%(size-1)!= 1 )) {
                        ia[j] = i+1;
                        ja[j] = i-(size-1);
                        aar[j] = stencilr[3];
                        j++;
                    }
                    
                    /* Down value */
                    if ( i >= size-1 ) {
                        ia[j] = i+1;
                        ja[j] = i+1-(size-1);
                        aar[j] = stencilr[4];
                        j++;
                    }
                    
                    /* SE value */
                    if (( i >= size-1 ) && ( (i+1)%(size-1)!= 0 )) {
                        ia[j] = i+1;
                        ja[j] = i+2-(size-1);
                        aar[j] = stencilr[5];
                        j++;
                    }
                    
                    /* NW value */
                    if (( i < (size-1 )*(size-2) ) && ( (i+1)%(size-1)!= 1 )){
                        ia[j] = i+1;
                        ja[j] = i + (size-1);
                        aar[j] = stencilr[6];
                        j++;
                    }
                    
                    
                    /* UP value */
                    if ( i < (size-1 )*(size-2) ) {
                        ia[j] = i+1;
                        ja[j] = i+1 + (size-1);
                        aar[j] = stencilr[7];
                        j++;
                    }
                    
                    /* NE value */
                    if (( i < (size-1 )*(size-2) ) && ( (i+1)%(size-1)!= 0 )){
                        ia[j] = i+1;
                        ja[j] = i+2+ (size-1);
                        aar[j] = stencilr[8];
                        j++;
                    }
                }
                break;
            } /* end case(9) */
            
            case ( 7 ):  {
                int sizeall = 7*(size-3)*(size-3)*(size-3)+36*(size-3)*(size-3)+60*(size-3)+32;
                plhs[0] = mxCreateDoubleMatrix(sizeall,1,mxREAL);
                plhs[1] = mxCreateDoubleMatrix(sizeall,1,mxREAL);
                plhs[2] = mxCreateDoubleMatrix(sizeall,1,mxREAL);
                
                ia = mxGetPr(plhs[0]);
                ja = mxGetPr(plhs[1]);
                aar = mxGetPr(plhs[2]);
                
                j = 0;
                for ( i = 0; i < (size -1)*(size -1)*(size -1); i++ ){
                    /* Center value */
                    ia[j] = i+1;
                    ja[j] = i+1;
                    aar[j] = stencilr[0];
                    j++;
                    
                    /* West value */
                    if ( i %(size-1)!= 0 ) {
                        ia[j] = i+1;
                        ja[j] = i;
                        aar[j] = stencilr[1];
                        j++;
                    }
                    
                    /* East value */
                    if ( (i+1)%(size-1)!= 0 ) {
                        ia[j] = i+1;
                        ja[j] = i+2;
                        aar[j] = stencilr[2];
                        j++;
                    }
                    
                    
                    /* South value */
                    if ( i % ((size-1)*(size-1)) >= size-1 ) {
                        ia[j] = i+1;
                        ja[j] = i+1-(size-1);
                        aar[j] = stencilr[3];
                        j++;
                    }
                    
                    
                    /* North value */
                    if ( i % ((size-1)*(size-1)) < (size-1)*(size-2) ) {
                        ia[j] = i+1;
                        ja[j] = i+1 + (size-1);
                        aar[j] = stencilr[4];
                        j++;
                    }
                    
                    /* Down value */
                    if ( i >= (size-1 )*(size-1) ) {
                        ia[j] = i+1;
                        ja[j] = i+1- (size-1)*(size-1);
                        aar[j] = stencilr[5];
                        j++;
                    }
                    
                    /* Up value */
                    if ( i < (size-1 )*(size-1)*(size-2) ) {
                        ia[j] = i+1;
                        ja[j] = i+1+ (size-1)*(size-1);
                        aar[j] = stencilr[6];
                        j++;
                    }
                }
                break;
            } /* end case(7) */
        }
    }
    else { /* Complex stencil */
        switch ( m )
        {
            case ( 3 ):  {
                plhs[0] = mxCreateDoubleMatrix(3*size-5,1,mxREAL);
                plhs[1] = mxCreateDoubleMatrix(3*size-5,1,mxREAL);
                plhs[2] = mxCreateDoubleMatrix(3*size-5,1,mxCOMPLEX);
                
                ia = mxGetPr(plhs[0]);
                ja = mxGetPr(plhs[1]);
                aar = mxGetPr(plhs[2]);
                aai = mxGetPi(plhs[2]);
                
                j = 0;
                for ( i = 0; i < size -1; i++ ){
                    ia[j] = i+1;
                    ja[j] = i+1;
                    aar[j] = stencilr[0];
                    aai[j] = stencili[0];
                    j++;
                    
                    if ( i - 1 >= 0 ) {
                        ia[j] = i+1;
                        ja[j] = i;
                        aar[j] = stencilr[1];
                        aai[j] = stencili[1];
                        j++;
                    }
                    
                    if ( i < size - 2 ) {
                        ia[j] = i+1;
                        ja[j] = i+2;
                        aar[j] = stencilr[2];
                        aai[j] = stencili[2];
                        j++;
                    }
                }
                break;
            }
            
            case ( 5 ):  {
                plhs[0] = mxCreateDoubleMatrix(5*size*size-14*size+9,1,mxREAL);
                plhs[1] = mxCreateDoubleMatrix(5*size*size-14*size+9,1,mxREAL);
                plhs[2] = mxCreateDoubleMatrix(5*size*size-14*size+9,1,mxCOMPLEX);
                
                ia = mxGetPr(plhs[0]);
                ja = mxGetPr(plhs[1]);
                aar = mxGetPr(plhs[2]);
                aai = mxGetPi(plhs[2]);
                
                j = 0;
                for ( i = 0; i < (size -1)*(size -1); i++ ){
                    /* Center value */
                    ia[j] = i+1;
                    ja[j] = i+1;
                    aar[j] = stencilr[0];
                    aai[j] = stencili[0];
                    j++;
                    
                    /* Left value */
                    if ( (i+1)%(size-1)!= 1 ) {
                        ia[j] = i+1;
                        ja[j] = i;
                        aar[j] = stencilr[1];
                        aai[j] = stencili[1];
                        j++;
                    }
                    
                    /* Right value */
                    if ( (i+1)%(size-1)!= 0 ) {
                        ia[j] = i+1;
                        ja[j] = i+2;
                        aar[j] = stencilr[2];
                        aai[j] = stencili[2];
                        j++;
                    }
                    
                    /* Down value */
                    if ( i >= size-1 ) {
                        ia[j] = i+1;
                        ja[j] = i+1-(size-1);
                        aar[j] = stencilr[3];
                        aai[j] = stencili[3];
                        j++;
                    }
                    
                    /* Up value */
                    if ( i < (size-1 )*(size-2) ) {
                        ia[j] = i+1;
                        ja[j] = i+1+ (size-1);
                        aar[j] = stencilr[4];
                        aai[j] = stencili[4];
                        j++;
                    }
                }
                break;
            } /* end case (5) */
            
            case ( 9 ):  {
                plhs[0] = mxCreateDoubleMatrix(3*(3*size-1)*(size-3)+16,1,mxREAL);
                plhs[1] = mxCreateDoubleMatrix(3*(3*size-1)*(size-3)+16,1,mxREAL);
                plhs[2] = mxCreateDoubleMatrix(3*(3*size-1)*(size-3)+16,1,mxCOMPLEX);
                
                ia = mxGetPr(plhs[0]);
                ja = mxGetPr(plhs[1]);
                aar = mxGetPr(plhs[2]);
                aai = mxGetPi(plhs[2]);
                aai = mxGetPi(plhs[2]);
                
                j = 0;
                for ( i = 0; i < (size -1)*(size -1); i++ ){
                    /* Center value */
                    ia[j] = i+1;
                    ja[j] = i+1;
                    aar[j] = stencilr[0];
                    aai[j] = stencili[0];
                    j++;
                    
                    /* Left value */
                    if ( (i+1)%(size-1)!= 1 ) {
                        ia[j] = i+1;
                        ja[j] = i;
                        aar[j] = stencilr[1];
                        aai[j] = stencili[1];
                        j++;
                    }
                    
                    /* Right value */
                    if ( (i+1)%(size-1)!= 0 ) {
                        ia[j] = i+1;
                        ja[j] = i+2;
                        aar[j] = stencilr[2];
                        aai[j] = stencili[2];
                        j++;
                    }
                    
                    /* SW value */
                    if (( i >= size-1 ) && ( (i+1)%(size-1)!= 1 )) {
                        ia[j] = i+1;
                        ja[j] = i-(size-1);
                        aar[j] = stencilr[3];
                        aai[j] = stencili[3];
                        j++;
                    }
                    
                    /* Down value */
                    if ( i >= size-1 ) {
                        ia[j] = i+1;
                        ja[j] = i+1-(size-1);
                        aar[j] = stencilr[4];
                        aai[j] = stencili[4];
                        j++;
                    }
                    
                    /* SE value */
                    if (( i >= size-1 ) && ( (i+1)%(size-1)!= 0 )) {
                        ia[j] = i+1;
                        ja[j] = i+2-(size-1);
                        aar[j] = stencilr[5];
                        aai[j] = stencili[5];
                        j++;
                    }
                    
                    /* NW value */
                    if (( i < (size-1 )*(size-2) ) && ( (i+1)%(size-1)!= 1 )){
                        ia[j] = i+1;
                        ja[j] = i + (size-1);
                        aar[j] = stencilr[6];
                        aai[j] = stencili[6];
                        j++;
                    }
                    
                    
                    /* UP value */
                    if ( i < (size-1 )*(size-2) ) {
                        ia[j] = i+1;
                        ja[j] = i+1 + (size-1);
                        aar[j] = stencilr[7];
                        aai[j] = stencili[7];
                        j++;
                    }
                    
                    /* NE value */
                    if (( i < (size-1 )*(size-2) ) && ( (i+1)%(size-1)!= 0 )){
                        ia[j] = i+1;
                        ja[j] = i+2+ (size-1);
                        aar[j] = stencilr[8];
                        aai[j] = stencili[8];
                        j++;
                    }
                }
                break;
            } /* end case(9) */
            
            case ( 7 ):  {
                int sizeall = 7*(size-3)*(size-3)*(size-3)+36*(size-3)*(size-3)+60*(size-3)+32;
                plhs[0] = mxCreateDoubleMatrix(sizeall,1,mxREAL);
                plhs[1] = mxCreateDoubleMatrix(sizeall,1,mxREAL);
                plhs[2] = mxCreateDoubleMatrix(sizeall,1,mxCOMPLEX);
                
                ia = mxGetPr(plhs[0]);
                ja = mxGetPr(plhs[1]);
                aar = mxGetPr(plhs[2]);
                aai = mxGetPi(plhs[2]);
                
                j = 0;
                for ( i = 0; i < (size -1)*(size -1)*(size -1); i++ ){
                    /* Center value */
                    ia[j] = i+1;
                    ja[j] = i+1;
                    aar[j] = stencilr[0];
                    aai[j] = stencili[0];
                    j++;
                    
                    /* West value */
                    if ( i %(size-1)!= 0 ) {
                        ia[j] = i+1;
                        ja[j] = i;
                        aar[j] = stencilr[1];
                        aai[j] = stencili[1];
                        j++;
                    }
                    
                    /* East value */
                    if ( (i+1)%(size-1)!= 0 ) {
                        ia[j] = i+1;
                        ja[j] = i+2;
                        aar[j] = stencilr[2];
                        aai[j] = stencili[2];
                        j++;
                    }
                    
                    
                    /* South value */
                    if ( i % ((size-1)*(size-1)) >= size-1 ) {
                        ia[j] = i+1;
                        ja[j] = i+1-(size-1);
                        aar[j] = stencilr[3];
                        aai[j] = stencili[3];
                        j++;
                    }
                    
                    
                    /* North value */
                    if ( i % ((size-1)*(size-1)) < (size-1)*(size-2) ) {
                        ia[j] = i+1;
                        ja[j] = i+1 + (size-1);
                        aar[j] = stencilr[4];
                        aai[j] = stencili[4];
                        j++;
                    }
                    
                    /* Down value */
                    if ( i >= (size-1 )*(size-1) ) {
                        ia[j] = i+1;
                        ja[j] = i+1- (size-1)*(size-1);
                        aar[j] = stencilr[5];
                        aai[j] = stencili[5];
                        j++;
                    }
                    
                    /* Up value */
                    if ( i < (size-1 )*(size-1)*(size-2) ) {
                        ia[j] = i+1;
                        ja[j] = i+1+ (size-1)*(size-1);
                        aar[j] = stencilr[6];
                        aai[j] = stencili[6];
                        j++;
                    }
                }
                break;
            } /* end case(7) */
        }
    }
    
    return;
}

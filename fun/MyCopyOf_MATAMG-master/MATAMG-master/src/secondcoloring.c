/*=========================================================
 * second_coloring.c
 * Second coloring scheme for AMG coarsening
 *
 *
 * This is a MEX-file for MATAMG.
 *=======================================================*/

#include <math.h>
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    
    double *C;
       
    int m = mxGetN(prhs[0]);
   
    double *ia, *ja, *aa,*is, *js;
    int is2ndFF = 0;
    int isNot2C = 0;
    int temp_index = 0;
    
    int i, j, k,jj;
    double a_ki;
    
    if(nrhs!=6) {
        mexErrMsgIdAndTxt("MATAMG:second_coloring:nrhs", "SECOND_COLORING(C,IA,JA,AA,IS,JS)");
    }
    
    /* get pointer  */
    C = mxGetPr(prhs[0]); 
    
    (void) plhs;
    
    ia = mxGetPr(prhs[1]);
    ja = mxGetPr(prhs[2]);
    aa = mxGetPr(prhs[3]);
    is = mxGetPr(prhs[4]);
    js = mxGetPr(prhs[5]);
        
    /* 
     * TODO make jsj
     * TODO compute Ci in advance 
     */
    
    for ( i = 1; i <= m; i++ ){
        if ( C[ i - 1 ] == 0.0) {    /* For an F-point i */
            is2ndFF = 0;
            /* Points that strongly influence to i */
            for ( j = (int)is[ i - 1 ]; j < (int)is[ i ]; j++ ) { 
               if ( C[ (int)js[ j - 1 ] - 1 ] == 0.0 && i != (int)js[ j - 1 ] ) {
                    isNot2C = 1;
                    for ( k = (int)is[ (int)js[ j - 1 ] - 1]; k < (int)is[ (int)js[ j - 1] ]; k++ ) {
                         if ( C[ (int)js[ k - 1 ] - 1] == 1.0 && (int)js[ k - 1 ] != i ) {
                            a_ki = 0.0;
                            for ( jj = (int)ia[ (int)js[ k - 1 ] - 1 ]; jj < (int)ia[ (int)js[ k - 1 ] ]; jj ++ ) {
                                if ( (int)ja[ jj - 1 ] == i ){
                                    a_ki = aa[ jj - 1 ];
                                }
                            }
                            if ( a_ki != 0.0 ){
                                isNot2C = 0;
                                break;
                            }   
                        }                         
                    } /* end loop k */
               
                    if ( isNot2C == 1 && is2ndFF == 1 ) { 
                        C[ i - 1 ] = 1.0;
                        C[ temp_index - 1 ] = 0.0;
                        
                        break;
                    }
                    
                    if ( isNot2C ==1 && is2ndFF ==0) {
                        is2ndFF = 1;
                        temp_index = (int)js[ j - 1 ];
                        C[ temp_index - 1 ]  = 1;
                    }
                }
            }
        }
    }
}

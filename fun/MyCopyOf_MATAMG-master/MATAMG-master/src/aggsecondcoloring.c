/*=========================================================
 * aggsecond_coloring.c
 * Second coloring scheme for aggressive coarsening

 * This is a MEX-file for MATAMG.
 *=======================================================*/

#include <math.h>
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    
    double *C;
       
    int m = mxGetN(prhs[0]);
   
    double *ia,*ja,*inbd, *jnbd, *is, *js;
    int is2ndFF = 0;
    int isNot2C = 0;
    int temp_index = 0;
    
    int i, j, k,jj;
    int jaj,jak,jnbdj, jnbdk;
    
    int *Ci, Ci_index;
    
    if(nrhs!=7) {
        mexErrMsgIdAndTxt("MATAMG:second_coloring:nrhs", "SECOND_COLORING(C,IA,JA,AA,IS,JS)");
    }
    
    /* get pointer  */
    C = mxGetPr(prhs[0]); 
    
    (void) plhs;
    
    ia = mxGetPr(prhs[1]); 
    ja = mxGetPr(prhs[2]);
    inbd = mxGetPr(prhs[3]);
    jnbd = mxGetPr(prhs[4]);
    is = mxGetPr(prhs[5]);
    js = mxGetPr(prhs[6]);
        
    for ( i = 1; i <= m; i++ ){
        if ( C[ i - 1 ] == 0.0) {    /* For an F-point i */
            Ci = mxCalloc((int)inbd[ i ] - (int)inbd[ i - 1 ] + (int)ia[ i ] - (int)ia[ i - 1 ],sizeof(int));
            
            Ci_index = 1;
            for ( j = (int)ia[ i - 1 ]; j < (int)ia[ i ]; j++ ){
                jaj = (int)ja[ j - 1 ];
                if (C[ jaj - 1 ] == 1 ){
                    Ci[ Ci_index - 1 ] = jaj;
                    Ci_index++;
                }
            }
            for ( j = (int) inbd[ i - 1 ]; j < (int)inbd[ i ]; j++ ){
                jnbdj = (int)jnbd[ j - 1 ];
                if (C[ jnbdj - 1 ] == 1 ){
                    Ci[ Ci_index - 1 ] = jnbdj;
                    Ci_index++;
                }
            }
                
            is2ndFF = 0;
            /* Points that strongly influence to i */
            for ( j = (int)ia[ i - 1 ]; j < (int)ia[ i ]; j++ ) { 
                jaj = (int)ja[ j - 1 ];
                if ( C[ (int)jaj - 1 ] == 0.0 && i != jaj ) {
                    isNot2C = 1;
                    for ( k = (int)ia[ jaj - 1]; k < (int)ia[ jaj ]; k++ ) {
                        jak = (int)ja[ k - 1 ];
                        if ( jak != i && jak != jaj ){
                            for ( jj = 1; jj < Ci_index; jj++){
                                if ( Ci[ jj - 1] == jak ){
                                    isNot2C = 0;
                                    break;
                                }
                            }
                        }
                        if (!isNot2C) break;
                    } /* end loop k */
                    if (isNot2C){
                        for ( k = (int)inbd[ jaj - 1]; k < (int)inbd[ jaj ]; k++ ) {
                            jnbdk = (int)jnbd[ k - 1 ];
                            for ( jj = 1; jj < Ci_index; jj++){
                                if ( Ci[ jj -1 ] == jnbdk ){
                                    isNot2C = 0;
                                    break;
                                }
                            }
                            if (!isNot2C) break;
                        } /* end loop k */
                    }
               
                    if ( isNot2C == 1 && is2ndFF == 1 ) { 
                        C[ i - 1 ] = 1.0;
                        C[ temp_index - 1 ] = 0.0;
                        break;
                    }
                    
                    if ( isNot2C ==1 && is2ndFF ==0) {
                        is2ndFF = 1;
                        temp_index = jaj;
                        C[ temp_index - 1 ]  = 1;
                    }
                }
            }
            mxFree(Ci);
        }
    }
}

/*=========================================================
 * longrangeamgp_mex.c
 *
 * This is a MEX-file for MATAMG.
 *=======================================================*/

#include <math.h>
#include "mex.h"

double getA(double *ia, double *ja, double *aa, int row, int col) {
    double a_i_j = 0.0;
    
    int j;
    for ( j = (int)ia[ row - 1 ]; j < (int)ia[ row ]; j++ ) {
        if ( (int)ja[ j - 1 ] == col ) {
            a_i_j = aa[ j - 1 ];
            break;
        }
    }        
    return a_i_j;    
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    double *C, *ia, *ja, *aa, *ip,*jp,*ap,*inbd,*jnbd, *carray;
    
    int pindex = 1;
    
    int m; /* Number of row of matrix A */
    
    int i,j,k,l;
    int *Ci, Ci_index;
    int jaj, jnbdj, jak, alreadyexist;
       
    if(nrhs!=10) {
        mexErrMsgIdAndTxt("MATAMG:longrangeamgp_mex:nrhs", "LONGRANGEAMGP(C,IA,JA,AA,IP,JP,AP,INBD,JNBD,CARRAY");
    }
       
    /* get pointer  */
    C = mxGetPr(prhs[0]);    
    
    ia = mxGetPr(prhs[1]);
    ja = mxGetPr(prhs[2]);
    aa = mxGetPr(prhs[3]);
    
    ip = mxGetPr(prhs[4]);
    jp = mxGetPr(prhs[5]);
    ap = mxGetPr(prhs[6]);
    
    inbd = mxGetPr(prhs[7]);
    jnbd = mxGetPr(prhs[8]);
    
    carray = mxGetPr(prhs[9]);
    
    m = mxGetN(prhs[0]);
    
    for ( i = 1; i <= m; i++ ){
        /*
         * mexPrintf("\ni = %d",i);
         */
        if ( C[ i - 1 ] == 1.0 ) {
            /*
             * mexPrintf(" C point.\n");
             */
            /*mexPrintf("[%d] is C point\n",i);*/
            ip[ pindex - 1 ] = i;
            jp[ pindex - 1 ] = carray[ i - 1 ];
            ap[ pindex - 1 ] = 1.0;
            pindex += 1;            
        }
        else {
            /*
             * mexPrintf(" F point.\n");
             */
            /*mexPrintf("[%d] is F point\n",i);*/
            Ci = mxCalloc((int)inbd[ i ] - (int)inbd[ i - 1 ] + (int)ia[ i ] - (int)ia[ i - 1 ],sizeof(int));
            
            Ci_index = 1;
                      
            /* Discrimate Ci,Dsi and Dwi among neighborhood points */
            /* 
             * TODO make jaj
             */
            for ( j = (int)ia[ i - 1 ]; j < (int)ia[ i ]; j++) {
                jaj = (int)ja[ j - 1 ];
                /*mexPrintf("nbd of %d  = %d\n",i,(int)ja[ j - 1 ]);*/
                if ( C[ jaj - 1 ] == 1.0 ) {
                    Ci[ Ci_index - 1 ] = jaj;
                    Ci_index += 1;
                    /*
                     * mexPrintf("\t[%d] => C_%d\n",jaj,i);
                     */
                }
            } /* end loop j */
            
            /* Find Ci among length 2 nbds */
            for ( j = (int)inbd[ i - 1 ]; j < (int)inbd[ i ]; j++){
                jnbdj = (int)jnbd[ j - 1 ];
                if ( C[ jnbdj - 1 ] == 1.0 ) {
                    Ci[ Ci_index - 1 ] = jnbdj;
                    Ci_index += 1;
                    /*
                     * mexPrintf("\t[%d] => C_%d\n",jnbdj,i);
                     */
                }
            }
            
            /* If no C point is found, check length 3 */
            if (Ci_index==1){
                for ( j = (int)inbd[ i - 1 ]; j < (int)inbd[ i ]; j++){
                    jnbdj = (int)jnbd[ j - 1 ];
                    for ( k = (int)ia[ jnbdj - 1 ] ; k < (int)ia[ jnbdj ];k++){
                        jak = (int)ja[ k - 1 ];
                        if ( C[ jak - 1 ] == 1.0 && jak != i) {
                            alreadyexist=0;
                            for (l = 1; l < Ci_index; l++){
                                if (jak == Ci[ l - 1 ]){
                                    alreadyexist=1;
                                    break;
                                }
                            }
                            if (!alreadyexist){
                                Ci[ Ci_index - 1 ] = jak;
                                Ci_index += 1;
                                /*
                                 * mexPrintf("\t[%d] => C_%d\n",jnbdj,i);
                                 */
                            }
                        }
                    }
                    
                }
            }
            
            if (Ci_index == 1)
                mexPrintf("no Ci points of %d\n",i);
            
            /*
             * mexPrintf("|Ci| : %d\n",Ci_index-1);
             */
            
            for ( j = 1; j < Ci_index; j++) {
                ip[ pindex - 1] = i;
                jp[ pindex - 1] = (int)carray[ Ci[ j - 1 ] - 1];
                ap[ pindex - 1] = 1.0/(Ci_index-1);         
                /*
                 * mexPrintf("Fine %d = Coarse %d\n",(int)Ci[ j - 1 ],(int)carray[ Ci[ j - 1 ] - 1]);
                 * mexPrintf("P(%d,%d) = %f\n",i,(int)carray[ Ci[ j - 1 ] - 1],  - 1.0/(Ci_index-1));                
                 */
                pindex += 1;
            }
            
            mxFree(Ci);            
        }
    }
    plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
    *mxGetPr(plhs[0])=(double)pindex;
    return;
}

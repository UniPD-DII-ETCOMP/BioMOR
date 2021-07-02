/*=========================================================
 * amgp_mex.c
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

double getS(double *is, double *js, int row, int col) {
    double a_i_j = 0.0;
    
    int j;
    for ( j = (int)is[ row - 1 ]; j < (int)is[ row ]; j++ ) {
        if ( (int)js[ j - 1 ] == col ) {
            a_i_j = 1.0;
            break;
        }
    }        
    return a_i_j;    
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    double *C, *ia, *ja, *aa, *ip,*jp,*ap,*is,*js, *carray;
    
    int pindex = 1;
    
    int m; /* Number of row of matrix A */
    
    int i,j,k, jaj,nnzrow;
    int *Ci, *Dsi, *Dwi, Ci_index, Dsi_index, Dwi_index;
    double a_i_k, a_k_j, *temp_Ci, bottom, bottom2;
       
    if(nrhs!=10) {
        mexErrMsgIdAndTxt("MATAMG:amgp_mex:nrhs", "AMGP_MEX(C,IA,JA,AA,IP,JP,AP,IS,JS,CARRAY");
    }
       
    /* get pointer  */
    C = mxGetPr(prhs[0]);    
    
    ia = mxGetPr(prhs[1]);
    ja = mxGetPr(prhs[2]);
    aa = mxGetPr(prhs[3]);
    
    ip = mxGetPr(prhs[4]);
    jp = mxGetPr(prhs[5]);
    ap = mxGetPr(prhs[6]);
    
    is = mxGetPr(prhs[7]);
    js = mxGetPr(prhs[8]);
    
    carray = mxGetPr(prhs[9]);
    
    m = mxGetN(prhs[0]);

    for ( i = 1; i <= m; i++ ){
        if ( C[ i - 1 ] == 1 ) {
            /*mexPrintf("[%d] is C point\n",i);*/
            ip[ pindex - 1 ] = i;
            jp[ pindex - 1 ] = carray[ i - 1 ];
            ap[ pindex - 1 ] = 1.0;
            pindex += 1;            
        }
        else {
            /*mexPrintf("[%d] is F point\n",i);*/
            nnzrow = (int)ia[ i ] - (int)ia[ i - 1 ];
            Ci = mxCalloc(nnzrow,sizeof(int));
            temp_Ci = mxCalloc(nnzrow,sizeof(double));
            Dsi = mxCalloc(nnzrow,sizeof(int));
            Dwi = mxCalloc(nnzrow,sizeof(int));
            
            Ci_index = 1;
            Dsi_index = 1;
            Dwi_index = 1;
                      
            bottom2 = 0;
            /* Discrimate Ci,Dsi and Dwi among neighborhood points */
            for ( j = (int)ia[ i - 1 ]; j < (int)ia[ i ]; j++) {
                jaj = (int)ja[ j - 1 ];
                /*mexPrintf("nbd of %d  = %d\n",i,(int)ja[ j - 1 ]);*/
                if ( C[ jaj - 1 ] == 1.0 ) {
                    Ci[ Ci_index - 1 ] = jaj;
                    ip[ pindex + Ci_index - 2 ] = i;
                    jp[ pindex + Ci_index - 2 ] = carray[ jaj - 1 ];
                    ap[ pindex + Ci_index - 2 ] = aa[ j - 1 ];
                    Ci_index += 1;
                    /*mexPrintf("\t[%d] => C_%d\n",(int)ja[ j - 1 ],i);*/
                }
                else {
                    if ( getS(is,js,i,jaj ) ){
                        Dsi[ Dsi_index - 1 ] = jaj;
                        Dsi_index += 1;
                        /*mexPrintf("\t[%d] => Ds_%d\n", (int)ja[ j - 1 ], i);*/
                    }
                    else if ( i != jaj ) {
                        Dwi[ Dwi_index - 1 ] =jaj;
                        bottom2 += aa[ j - 1 ];
                        Dwi_index += 1;
                        /*mexPrintf("\t[%d] => Dw_%d\n",(int)ja[ j - 1 ],i);*/
                    }
                    else {
                        bottom2 += aa[ j - 1 ];
                        /*mexPrintf("\t[%d] = %d\n",(int)ja[ j - 1 ],i);*/
                    }                    
                }
            } /* end loop j */
            
            /* calculate weight P(i,j)*/
            for ( k = 1; k < Dsi_index; k++ ) {
                bottom = 0.0;
                a_i_k = getA(ia,ja,aa,i,Dsi[ k - 1 ]);
                for ( j = 1; j < Ci_index; j++ ) {
                    a_k_j = getA(ia,ja,aa,Ci[ j - 1 ], Dsi[ k - 1 ] );
                    if ( a_k_j != 0.0 ) {
                        temp_Ci[ j - 1 ] += a_i_k*a_k_j;
                        bottom += a_k_j;
                    }
                }
                
                for ( j = 1; j < Ci_index; j++ ) {
                     ap[ pindex + j - 2 ] += temp_Ci[ j - 1 ]/bottom;
                     temp_Ci[ j - 1 ] = 0;
                }
            }
            
            for ( j = 0; j < Ci_index; j++) {
                ap[ pindex + j - 1] = - ap[ pindex + j -1 ]/bottom2;                
            }
            
            pindex += Ci_index - 1;
            mxFree(Ci);
            mxFree(Dsi); 
            mxFree(Dwi);
            mxFree(temp_Ci);
        }
    }
    plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
    *mxGetPr(plhs[0])=(double)pindex;
    return;
}

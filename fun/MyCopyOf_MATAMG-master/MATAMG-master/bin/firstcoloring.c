/*=========================================================
 * first_coloring.c
 * First coloring scheme for AMG coarsening
 *
 *  nofc        -1  : assigned
 *              o.w : Number of connections
 *
 *  prev_index  -1  : deleted or unset
 *               0  : the first node
 *              o.w.: index to the previous node
 *
 *  next_index  -1  : the last node
 *               0  : deleted or unset
 *              o.w.: index to the next node
 *
 *  first_index  0  : no node
 *              o.w : index to the first node
 * 
 *
 * This is a MEX-file for MATAMG.
 *=======================================================*/

#include <math.h>
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double *C;
    int m;
    int *prev_index,*next_index, *first_index;
    double *it, *jt, *is, *js, *nofc;
    int max_lambda = 0;
    int lambda, lambdam1;
    
    int i;
    int assigned_index = 0;     /* Number of assigned nodes*/
    int temp_index;
    int max_index;
    int l, j, k; 
    int jtjm1,jskm1;
    
    if(nrhs!=6) {
        mexErrMsgIdAndTxt("MATAMG:first_coloring:nrhs", "FIRST_COLORING(C,IT,JT,IS,JS)");
    }
    
    (void) plhs;
    
    /* get pointer  */
    C = mxGetPr(prhs[0]);    
     
    m = mxGetN(prhs[0]);
   
    prev_index = mxCalloc(m,sizeof(int));
    next_index = mxCalloc(m,sizeof(int));
    first_index = mxCalloc(m,sizeof(int));
    
    
    it = mxGetPr(prhs[1]);
    jt = mxGetPr(prhs[2]);
    is = mxGetPr(prhs[3]);
    js = mxGetPr(prhs[4]);
    nofc = mxGetPr(prhs[5]);
    
    /* Initialize first_index */
    for ( i = 0; i < m; i++ ){ 
        first_index[ i ] = 0;
        prev_index[ i ] = -1;
        next_index[ i ] = 0;
    }
    
    /* Initialize nofc and double linked list, and update max_lambda */
    for ( i = 0; i < m; i++ ) {
        lambda = (int)nofc[ i ];
        lambdam1 = lambda - 1;
        /*
         * TODO first_index[ lambda - 1 ] x 3 lambda -1 x 5
         */
        if (lambda > 0) {
            if (max_lambda < lambda) {    /* Update maximum lambda value */
                max_lambda = lambda;
            }
            
            if ( first_index[ lambdam1 ] == 0 ) {   /* First new node */
                first_index[ lambdam1 ] = i + 1;
                prev_index[ i ] = 0;
                next_index[ i ] = -1;                
            }
            else {                                  /* Add new node */
                prev_index[ first_index[ lambdam1 ] - 1 ] = i + 1;
                next_index[ i ] = first_index[ lambdam1 ];
                
                prev_index[ i ] = 0;
                first_index[ lambdam1 ] = i + 1;
            }
        }
    }
    
    
    /* First coloring */
    while ( assigned_index != m ) {
        if ( max_lambda > 0 ) {
            max_index = first_index[ max_lambda - 1 ];
            temp_index = max_index;
            
            /*
             * TODO next_index[ temp_index - 1 ] x 4
             */
            while ( next_index[ temp_index - 1 ]  != -1 ) {
                if ( max_index > next_index[ temp_index - 1 ] ) {
                    max_index = next_index[ temp_index - 1 ];
                }
                temp_index = next_index[ temp_index - 1 ];                
            } /* end while */                    
        } /* end if ( max_lambda > 0 ) */
        else {
            for ( l = 1; l <= m; l++ ) {
                if ( nofc[ l - 1 ] != -1 )      
                    C[ l - 1 ] = 1;                
            } /* end for loop l*/
            break;
        } /* end else */
        
        C[ max_index - 1 ] = 1;
        nofc[ max_index -1 ] = -1;
        assigned_index += 1;
        
        /* delete node from the list */
        if ( next_index[ max_index - 1 ] == -1 ) {     /* Last index of list */
            if ( prev_index[ max_index - 1 ] == 0 ) {    /* First index of list */
                first_index[ max_lambda - 1 ] = 0;
                while ( first_index[ max_lambda - 1] == 0 ) {
                    max_lambda -= 1;
                    if ( max_lambda == 0 )
                        break;
                } /* end while*/
            } /* end if ( prev...) */
            else {                                      /* Last index of list of size > 1 */
                next_index[ prev_index[ max_index - 1 ] -1 ] = -1;
            } /* end else*/
        } /* end if ( next... ) */
        else {
            if ( prev_index[ max_index - 1 ] == 0 ) {   /* First index of list of size > 1 */
                first_index[ max_lambda - 1 ] = next_index[ max_index - 1];
                prev_index[ next_index[ max_index - 1 ] - 1 ] = 0;
            } /* end if ( prev...) */
            else {      /* index in the middle of list of size > 1 */
                prev_index[ next_index[ max_index - 1 ] -1 ] = prev_index[ max_index - 1 ];
                next_index[ prev_index[ max_index - 1 ] -1 ] = next_index[ max_index - 1 ];
            } /* end else*/
        } /* end else */
        
        /* Make pointers of max_index Null*/
        prev_index[ max_index - 1 ] = -1;
        next_index[ max_index - 1 ] = 0;
        
        /* Set unassigned points which strongly depend on max_index to be F point */
        for ( j = (int)it[ max_index - 1 ]; j < (int)it[ max_index ]; j++ ) {
            jtjm1 = (int)jt[ j - 1 ] - 1;
            if ( nofc[ jtjm1 ] > 0 ) {
                /* Delete node jt(j)*/
                if ( next_index[ jtjm1 ] == -1 ) { /* Last index of list */
                    if ( prev_index[ jtjm1 ] == 0 ) { /* First index of list */
                        first_index[ (int)nofc[ jtjm1 ] - 1] = 0;
                        while ( first_index[ max_lambda - 1 ] == 0 ) {
                            max_lambda -= 1;
                            if ( max_lambda == 0 )
                                break;
                        } /* end while */
                    } /* end if */
                    else { /* Last index of list of size > 1 */
                        next_index[ prev_index[ jtjm1 ] - 1 ] = -1;
                    }
                } /* end if (next...) */
                else {
                    if (prev_index[ jtjm1 ] == 0 ) { /* First index of list of size > 1 */
                        first_index[ (int)nofc[ jtjm1 ] - 1] = next_index[ jtjm1 ];
                        prev_index[ next_index[ jtjm1 ] - 1] = 0;
                    }
                    else {
                        prev_index[ next_index[ jtjm1 ] - 1 ] = prev_index[ jtjm1 ];
                        next_index[ prev_index[ jtjm1 ] - 1 ] = next_index[ jtjm1 ];
                    }
                } /* end else */
                
                /* Make pointers of max index NULL */
                prev_index[ jtjm1 ] = -1;
                next_index[ jtjm1 ] =0;
                
                /* Assign jt(j) to be F point */
                nofc[ jtjm1 ] = - 1;
                assigned_index += 1;
                
                
                /* Check nbds of jt(j) */
                for ( k = (int)is[ (int)jt[ j - 1 ] - 1 ]; k < (int)is[ (int)jt[ j -1 ] ]; k++ ) {
                    jskm1 = (int)js[ k - 1 ] - 1;
                    if ( nofc[ jskm1 ] > 0 ) {
                        /* Delete js(k) from list of nofc(js(k))*/
                        if ( next_index[ jskm1 ] == -1){ /* Last index of list*/
                            if ( prev_index[ jskm1 ] == 0 ) { /* First index of list */
                                first_index[ (int)nofc[ jskm1 ] - 1 ] = 0;
                            }
                            else {
                                next_index[ prev_index[ jskm1 ] - 1 ] = -1;
                            }
                        }
                        else {
                            if ( prev_index[ jskm1 ] == 0 ) {
                                first_index[ (int)nofc[ jskm1 ] - 1 ] = next_index[ jskm1 ];
                                prev_index[ next_index[ jskm1 ] - 1 ] = 0;
                            }
                            else {
                                prev_index[ next_index[ jskm1 ] - 1 ] = prev_index[ jskm1 ];
                                next_index[ prev_index[ jskm1 ] - 1 ] = next_index[ jskm1 ];
                            }
                        }
                        
                        /* Increase the number of connection by 1 */
                        nofc[ jskm1 ] += 1;
                        
                        /* Update max_lambda */
                        if ( nofc[ jskm1 ] > max_lambda ){
                            max_lambda = nofc[ jskm1 ];
                        }
                        /* Insert index to new list */
                        if ( first_index[ (int)nofc[ jskm1 ] - 1 ] == 0 ) {
                            first_index[ (int)nofc[ jskm1 ] - 1 ] = jskm1 + 1;
                            prev_index[ jskm1 ] = 0;
                            next_index[ jskm1 ] = -1;
                        }
                        else {
                            prev_index[ first_index[ (int)nofc[ jskm1 ] - 1 ] - 1 ] = jskm1 + 1;
                            next_index[ jskm1 ] = first_index[ (int)nofc[ jskm1 ] - 1 ];
                            
                            prev_index[ jskm1 ] = 0;
                            first_index[ (int)nofc[ jskm1 ] - 1 ] = jskm1 + 1;
                        }
                    } /* end if (nofc...) */
                } /* end loop k */            
            } /* end if */
        } /* end loop j*/
    } /* end while */
    mxFree(prev_index);
    mxFree(next_index);
    mxFree(first_index);
}

/*=========================================================
 * aggcoarsenging_mex.c
 * aggresive coarsening for AMG coarsening
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
    
    
    double *C, *C2;
    int m;
    int *nofc,*prev_index,*next_index, *first_index;
    double *it, *jt, *is, *js;
    int max_lamda = 0;
    int lamda;
    int path;
    
    int i;
    int assigned_index = 0;     /* Number of assigned nodes*/
    int temp_index;
    int max_index;
    int l, j, k; 
    
    if(nrhs!= 6) {
        mexErrMsgIdAndTxt("MATAMG:aggcoarsening_mex:nrhs", "AGGCOARSENING_MEX(C2,C,IT,JT,IS,JS)");
    }
    
    (void) plhs;
    
    /* get pointer  */
    C = mxGetPr(prhs[0]);    
    C2 = mxGetPr(prhs[1]); 
    m = mxGetN(prhs[0]);
   
    nofc       = mxCalloc(m,sizeof(int));
    prev_index = mxCalloc(m,sizeof(int));
    next_index = mxCalloc(m,sizeof(int));
    first_index = mxCalloc(m,sizeof(int));
    
    
    it = mxGetPr(prhs[2]);
    jt = mxGetPr(prhs[3]);
    is = mxGetPr(prhs[4]);
    js = mxGetPr(prhs[5]);
         
    
    
    /* Initialize first_index */
    for ( i = 1; i <= m; i++ ){ 
   /* for ( i = m; i >= 1; i-- ){*/
        first_index[ i - 1 ] = 0;
        prev_index[ i - 1 ] = -1;
        next_index[ i - 1 ] = 0;
    }
    
    /* Initialize nofc and double linked list, and update max_lamda */
    for ( i = 1; i <= m; i++ ) {
    /*for ( i = m; i >= 1; i-- ){*/
        lamda = (int)it[ i ] - (int)it[ i - 1 ] + (int)is[ i ] - (int)is[ i - 1 ];
        
        /*
         * mexPrintf("lambda[%d] = %d\n",i,lamda);
        */
        
        nofc[ i - 1 ] = lamda;
        if (lamda > 0) {
            if (max_lamda < lamda) {    /* Update maximum lamda value */
                max_lamda = lamda;
            }
            
            if ( first_index[ lamda -1 ] == 0 ) {   /* First new node */
                first_index[ lamda -1 ] = i;
                prev_index[ i - 1 ] = 0;
                next_index[ i - 1 ] = -1;                
            }
            else {                                  /* Add new node */
                prev_index[ first_index[ lamda - 1 ] - 1 ] = i;
                next_index[ i - 1 ] = first_index[ lamda - 1 ];
                
                prev_index[ i - 1 ] = 0;
                first_index[ lamda - 1 ] = i;
            }
        }
    }
   
    
    
    
    /* First coloring */
    while ( assigned_index != m ) {
        if ( max_lamda > 0 ) {
            max_index = first_index[ max_lamda - 1 ];
            temp_index = max_index;
                    
            while ( next_index[ temp_index - 1 ]  != -1 ) {
                if ( max_index > next_index[ temp_index - 1 ] ) {
                    max_index = next_index[ temp_index - 1 ];
                }
                temp_index = next_index[ temp_index - 1 ];                
            } /* end while */                    
        } /* end if ( max_lamda > 0 ) */
        else {
            for ( l = 1; l <= m; l++ ) {
                if ( nofc[ l - 1 ] != -1 && C[ l - 1 ] == 1)      
                    C2[ l - 1 ] = 1;                
            } /* end for loop l*/
            break;
        } /* end else */
        
        /*
         * mexPrintf("%d => C\n",max_index);
        */
        
        C2[ max_index - 1 ] = 1;
        nofc[ max_index -1 ] = -1;
        assigned_index += 1;
        
        /* delete node from the list */
        if ( next_index[ max_index - 1 ] == -1 ) {     /* Last index of list */
            if ( prev_index[ max_index - 1 ] == 0 ) {    /* First index of list */
                first_index[ max_lamda - 1 ] = 0;
                while ( first_index[ max_lamda - 1] == 0 ) {
                    max_lamda -= 1;
                    if ( max_lamda == 0 )
                        break;
                } /* end while*/
            } /* end if ( prev...) */
            else {                                      /* Last index of list of size > 1 */
                next_index[ prev_index[ max_index - 1 ] -1 ] = -1;
            } /* end else*/
        } /* end if ( next... ) */
        else {
            if ( prev_index[ max_index - 1 ] == 0 ) {   /* First index of list of size > 1 */
                first_index[ max_lamda - 1 ] = next_index[ max_index - 1];
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
            if ( nofc[ (int)jt[ j - 1 ] - 1 ] > 0 ) {
                /*
                 * mexPrintf("Deleting %d\n",(int)jt[ j - 1 ]);
                */
                
                /* Delete node jt(j)*/
                if ( next_index[ (int)jt[ j - 1 ] - 1 ] == -1 ) { /* Last index of list */
                    if ( prev_index[ (int)jt[ j - 1 ] - 1 ] == 0 ) { /* First index of list */
                        first_index[ nofc[ (int)jt[ j - 1 ] - 1 ] - 1] = 0;
                        while ( first_index[ max_lamda - 1 ] == 0 ) {
                            max_lamda -= 1;
                            if ( max_lamda == 0 )
                                break;
                        } /* end while */
                    } /* end if */
                    else { /* Last index of list of size > 1 */
                        next_index[ prev_index[ (int)jt[ j - 1 ] - 1 ] - 1 ] = -1;
                    }
                } /* end if (next...) */
                else {
                    if (prev_index[ (int)jt[ j- 1 ] - 1 ] == 0 ) { /* First index of list of size > 1 */
                        first_index[ nofc[ (int)jt[ j - 1 ] - 1 ] - 1] = next_index[ (int)jt[ j - 1] - 1 ];
                        prev_index[ next_index[ (int)jt[ j - 1 ] - 1] - 1] = 0;
                    }
                    else {
                        prev_index[ next_index[ (int)jt[ j - 1 ] - 1 ] - 1 ] = prev_index[ (int)jt[ j - 1 ] - 1 ];
                        next_index[ prev_index[ (int)jt[ j - 1 ] - 1 ] - 1 ] = next_index[ (int)jt[ j - 1 ] - 1 ];
                    }
                } /* end else */
                
                /* Make pointers of max index NULL */
                prev_index[ (int)jt[ j - 1 ] - 1 ] = -1;
                next_index[ (int)jt[ j - 1 ] - 1 ] =0;
                
                /* Assign jt(j) to be F point */
                nofc[ (int)jt[ j - 1 ] - 1 ] = - 1;
                assigned_index += 1;
                
                
                /* Check nbds of jt(j) */
                for ( k = (int)is[ (int)jt[ j - 1 ] - 1 ]; k < (int)is[ (int)jt[ j -1 ] ]; k++ ) {
                    if ( nofc[ (int)js[ k - 1 ] - 1 ] > 0 ) {
                        /* Delete js(k) from list of nofc(js(k))*/
                        if ( next_index[ (int)js[ k - 1 ] - 1 ] == -1){ /* Last index of list*/
                            if ( prev_index[ (int)js[ k - 1 ] - 1 ] == 0 ) { /* First index of list */
                                first_index[ nofc[ (int)js[ k - 1 ] - 1 ] - 1 ] = 0;
                            }
                            else {
                                next_index[ prev_index[ (int)js[ k - 1 ] - 1 ] - 1 ] = -1;
                            }
                        }
                        else {
                            if ( prev_index[ (int)js[ k - 1 ] - 1 ] == 0 ) {
                                first_index[ nofc[ (int)js[ k - 1 ] - 1 ] - 1 ] = next_index[ (int)js[ k - 1 ] - 1 ];
                                prev_index[ next_index[ (int)js[ k - 1 ] - 1 ] - 1 ] = 0;
                            }
                            else {
                                prev_index[ next_index[ (int)js[ k - 1 ] - 1 ] - 1 ] = prev_index[ (int)js[ k - 1 ] - 1 ];
                                next_index[ prev_index[ (int)js[ k - 1 ] - 1 ] - 1 ] = next_index[ (int)js[ k - 1 ] - 1 ];
                            }
                        }
                        
                        /* Increase the number of connection by 1 */
                        nofc[ (int)js[ k - 1 ] - 1 ] += 1;
                        /*
                         * mexPrintf("%d ++ \n",(int)js[ k - 1 ]);
                         */
                        /* Update max_lamda */
                        if ( nofc[ (int)js[ k - 1 ] - 1 ] > max_lamda ){
                            max_lamda = nofc[ (int)js[ k - 1 ] - 1 ];
                        }
                        /* Insert index to new list */
                        if ( first_index[ nofc[ (int)js[ k - 1 ] - 1 ] - 1 ] == 0 ) {
                            first_index[ nofc[ (int)js[ k - 1 ] - 1 ] - 1 ] = (int)js[ k - 1 ];
                            prev_index[ (int)js[ k - 1 ] - 1 ] = 0;
                            next_index[ (int)js[ k - 1 ] - 1 ] = -1;
                        }
                        else {
                            prev_index[ first_index[ nofc[ (int)js[ k - 1 ] - 1 ] - 1 ] - 1 ] = (int)js[ k - 1 ];
                            next_index[ (int)js[ k - 1 ] - 1 ] = first_index[ nofc[ (int)js[ k - 1 ] - 1 ] - 1 ];
                            
                            prev_index[ (int)js[ k - 1 ] - 1 ] = 0;
                            first_index[ nofc[ (int)js[ k - 1 ] - 1 ] - 1 ] = (int)js[ k - 1 ];
                        }
                    } /* end if (nofc...) */
                } /* end loop k */            
            } /* end if */
        } /* end loop j*/
    } /* end while */
    
    mxFree(nofc);
    mxFree(prev_index);
    mxFree(next_index);
    mxFree(first_index);
}

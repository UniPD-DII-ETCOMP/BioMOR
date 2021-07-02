#ifndef matamg
#include "mex.h"
#include <math.h>
#endif

#define SQR(a) (a==0.0? 0.0 : a*a)
#define SIGN(a, b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define FMAX(a, b) ( a > b ? a : b )

/* Get real value of matrix */
double getR(const mxArray *A, int i, int j) {
    double *pr;
    size_t m;
    
    pr = mxGetPr(A);
    m = mxGetM(A);
    return pr[m*j+i];
}

/* Get imaginary value of matrix */
double getI(const mxArray *A, int i, int j) {
    double *pti;
    size_t m;
    pti = mxGetPi(A);
    m = mxGetM(A);
    return pti[m*j+i];
}

/* Set real value of matrix */
void setR(mxArray *A, int i, int j, double v) {
    double *pr;
    size_t m;
    pr = mxGetPr(A);
    m = mxGetM(A);
    pr[m*j+i] = v;
}

/* Set imaginary value of matrix */
void setI(mxArray *A, int i, int j, double v) {
    double *pti;
    size_t m;
    pti = mxGetPi(A);
    m = mxGetM(A);
    pti[m*j+i] = v;
}

/* multiplication for real valued matrices */
mxArray* times( const mxArray *A, const mxArray *B) {
    double *ar, *br, *cr;
    size_t m, n, l;
    int i, j, k;
    mxArray *C;
    
    if (mxGetN(A) != mxGetM(B) ) {
        mexErrMsgTxt("Inner matrix dimensions must agree.");
    }
    
    C = mxCreateDoubleMatrix((int)mxGetM(A), (int)mxGetN(B), mxREAL);
    
    ar = mxGetPr(A);
    br = mxGetPr(B);
    cr = mxGetPr(C);
    
    m = mxGetM(A);
    n = mxGetN(B);
    l = mxGetN(A);
    
    for( i = 0; i < m; i++ ){
        for( j = 0; j < n; j++ ) {
            for( k = 0; k < l; k++ ) {
                cr[ m*j+i ] += ar[ m*k+i ]*br[ m*j+k ];
            }
        }
    }
    return C;
}

/* A = A - B */
void minus( mxArray *A, const mxArray *B ) {
    
    double *apr, *bpr;
    int i;
    
    if (mxGetM(A) != mxGetM(B) || mxGetN(A) != mxGetN(B) ) {
        mexErrMsgTxt("Inner matrix dimensions must agree.");
    }
    
    
    apr = mxGetPr(A);
    bpr = mxGetPr(B);
    
    
    for ( i = 0; i < mxGetM(A)*mxGetN(A); i++ ){
        apr[i] -= bpr[i];
    }
}

/* A = A/numerator */
void divide( mxArray *A, double numerator ) {
    double *apr = mxGetPr(A);
    int i;
    
    for ( i = 0; i < mxGetM(A)*mxGetN(A); i++ ){
        apr[i] /= numerator;
    }
}
/* multiplication for real valued matrices */
void times2( const mxArray *A, const mxArray *B, mxArray *C) {
    
    double *ar, *br, *cr;
    size_t m, n, l;
    int i, j, k;
    
    if (mxGetN(A) != mxGetM(B) ) {
        mexErrMsgTxt("Inner matrix dimensions must agree.");
    }
    
    ar = mxGetPr(A);
    br = mxGetPr(B);
    cr = mxGetPr(C);
    
    m = mxGetM(A);
    n = mxGetN(B);
    l = mxGetN(A);
    
    for( i = 0; i < m; i++ ){
        for( j = 0; j < n; j++ ) {
            for( k = 0; k < l; k++ ) {
                cr[ m*j+i ] += ar[ m*k+i ]*br[ m*j+k ];
            }
        }
    }
}

/* build diagonal matrix of a reciprocal value. */
mxArray* idiag( const mxArray *d) {
    
    mxArray *iD;
    size_t m;
    double *idpr;
    double *dpr;
    int i;
    
    m = ( mxGetM(d) > mxGetN(d) ? mxGetM(d): mxGetN(d));
    if ( !( mxGetM(d)== 1||mxGetN(d) == 1 ) ) {
        mexErrMsgTxt("Input must be a vector.");
    }
    iD = mxCreateDoubleMatrix((int)m, (int)m, mxREAL);
    
    idpr = mxGetPr(iD);
    dpr = mxGetPr(d);
    
    for ( i = 0; i < m; i++ ) {
        idpr[ (m+1)*i ] = 1.0/dpr[ i ];
    }
    return iD;
}

/* return submatirx. */
mxArray* sub( const mxArray *A, int row1, int row2, int col1, int col2) {
    size_t m, n;
    double *Apr, *subApr;
    
    int i, j;
    mxArray *subA;
    
    m = mxGetM(A);
    n = mxGetN(A);
    
    if ( row1 > row2 || col1 > col2 || row1 < 0 || row2 > m-1 || col1 < 0 || col2 > n-1 ) {
        mexErrMsgTxt("Wrong indexes are used.");
    }
    
    
    subA = mxCreateDoubleMatrix(row2-row1+1, col2-col1+1, mxREAL);
    
    Apr = mxGetPr(A);
    subApr = mxGetPr(subA);
    
    for ( i = row1; i <= row2 ; i++ ) {
        for ( j = col1; j <= col2; j++ ) {
            subApr[ (row2-row1+1)*(j-col1) +(i-row1)] = Apr[ m*j+i];
        }
    }
    return subA;
}

/* return submatirx( only for rbamg ) . */
mxArray* subrbamg( const mxArray *A, const mxArray *R, double d, int row1, int row2, int col1, int col2) {
    size_t m, n;
    
    
    double *Apr, *Rpr, *subApr;
    int i, j;
    
    mxArray *subA;
    
    m = mxGetM(A);
    n = mxGetN(A);
    
    if ( row1 > row2 || col1 > col2 || row1 < 0 || row2 > m-1 || col1 < 0 || col2 > n-1 ) {
        mexErrMsgTxt("Wrong indexes are used.");
    }
    
    subA = mxCreateDoubleMatrix(row2-row1+1, col2-col1+1, mxREAL);
    
    Apr = mxGetPr(A);
    Rpr = mxGetPr(R);
    subApr = mxGetPr(subA);
    
    
    for ( i = row1; i <= row2 ; i++ ) {
        for ( j = col1; j <= col2; j++ ) {
            subApr[ (row2-row1+1)*(j-col1) +(i-row1)] = Apr[ m*j+i] -Rpr[ m*j+i ]/d;
        }
    }
    return subA;
}



/* return submatirx( only for ibamg ) which is used in Least squares process. */
mxArray* subibamg( const mxArray *A, double *a_i_k, int *Fi, int length) {
    size_t m, n;
    
    double *Apr, *subApr;
    
    int i, j;
    
    mxArray *subA;
    
    m = mxGetM(A);
    n = mxGetN(A);
    
    subA = mxCreateDoubleMatrix((int)m, 1, mxREAL);
    
    Apr = mxGetPr(A);
    subApr = mxGetPr(subA);
    
    for ( i = 0; i < m ; i++ ) {
        for ( j = 0;  j < length; j++ ) {
            subApr[ i ] += Apr[ m*(Fi[ j ] - 1) +i ] * a_i_k[j];
        }
    }
    return subA;
}

/* return submatirx( only for ibamg ) which is used in Least squares process. */
mxArray* subsibamg( const mxArray *A, int col) {
    size_t m, n;
    double *Apr, *subApr;
    mxArray *subA;
    int i;
    
    m = mxGetM(A);
    n = mxGetN(A);
    
    
    subA = mxCreateDoubleMatrix((int)m, 1, mxREAL);
    
    Apr = mxGetPr(A);
    subApr = mxGetPr(subA);
    
    for ( i = 0; i < m ; i++ ) {
        subApr[ i ] = Apr[m*(col - 1) + i];
        
    }
    return subA;
}

mxArray* subcol( const mxArray *A, int *index, int length ) {
    size_t m, n;
    double *Apr, *subApr;
    mxArray *subA;
    int i, j;
    
    m = mxGetM(A);
    n = mxGetN(A);
    
    
    subA = mxCreateDoubleMatrix((int)m, length, mxREAL);
    
    Apr = mxGetPr(A);
    subApr = mxGetPr(subA);
    
    
    for ( i = 0; i < m; i++ ) {
        for ( j = 0; j < length; j++ ) {
            subApr[ m*j+i ] = Apr[ m*(index[j]-1)+i];
        }
    }
    return subA;
    
}
/* Copy matrix A to B */
void copy( const mxArray *A, mxArray *B) {
    
    int i;
    double *ar, *br;
    if( mxGetM(A) != mxGetM(B) || mxGetN(A) != mxGetN(B) ) {
        mexErrMsgTxt("Size of matrix you want to copy is different from that of orginal matrix.");
    }
    ar = mxGetPr(A);
    br = mxGetPr(B);
    
    for ( i = 0; i < mxGetM(A)*mxGetN(A); i++){
        br[i] = ar[i];
    }
}



double pythag(double a, double b) {
    
    double absa = fabs(a);
    double absb = fabs(b);
    
    if(absa > absb) {
        return absa*sqrt(1.0+SQR(absb/absa));
    }
    else {
        return (absb == 0.0? 0.0:absb*sqrt(1.0+SQR(absa/absb)));
    }
}

mxArray* transpose( const mxArray *A) {
    
    double *ar, *atr;
    mxArray *AT;
    
    int i, j;
    
    ar = mxGetPr(A);
    
    AT = mxCreateDoubleMatrix((int)mxGetN(A), (int)mxGetM(A), mxREAL);
    atr = mxGetPr(AT);
    
    
    for ( i = 0; i < mxGetN(A); i++ ) {
        for ( j = 0; j < mxGetM(A); j++ ) {
            atr[ mxGetN(A)*j + i ] = ar[ mxGetM(A)*i + j ];
        }
    }
    
    return AT;
    
}



void svdcmp( const mxArray *A, mxArray *U, mxArray *V, mxArray *w) {
    
    int flag, i, its, j, jj, k, l, nm;
    double anorm, c, f, g, h, s, scale, x, y, z;
    
    size_t m = mxGetM(A);
    size_t n = mxGetN(A);
    double *wr = mxGetPr(w);
    
    
    double *rv1;
    copy(A, U);
    
    
    rv1 = mxCalloc(n, sizeof(double));
    g = scale = anorm = 0.0;
    
    for (i = 0; i < n; i++ ) {
        l = i + 1;
        rv1[ i ] = scale * g;
        g = s = scale = 0.0;
        if ( i < m ) {
            for ( k = i; k < m; k++ ) {
                scale += fabs(getR(U, k, i));
            } /* end loop k*/
            
            if ( scale != 0.0 ) {
                for ( k = i; k < m; k++ ) {
                    setR(U, k, i, getR(U, k, i)/scale);
                    s += getR(U, k, i)*getR(U, k, i);
                } /* end loop k */
                
                f = getR(U, i, i);
                g = -SIGN(sqrt(s), f);
                h = f*g-s;
                setR(U, i, i, f-g);
                
                for ( j = l; j < n; j++ ) {
                    for ( s = 0.0, k = i; k < m; k++ ) {
                        s += getR(U, k, i)*getR(U, k, j);
                    }
                    f = s/h;
                    for ( k = i; k < m; k++ ) {
                        setR(U, k, j, getR(U, k, j) + f * getR(U, k, i));
                    }
                } /* end loop j */
                
                for ( k = i; k < m; k++ ) {
                    setR(U, k, i, getR(U, k, i)*scale);
                }
            } /* end if ( scale !=... */
        } /* end if ( i < m ) */
        
        wr[i] = scale * g;
        g = s = scale = 0.0;
        
        if ( i < m && i != n - 1 ) {
            for ( k = l; k < n; k++ ) {
                scale += fabs(getR(U, i, k));
            } /* end loop k */
            
            if ( scale != 0 ) {
                for ( k = l; k < n; k++ ) {
                    setR(U, i, k, getR(U, i, k)/scale);
                    s+= getR(U, i, k) * getR(U, i, k);
                } /* end loop k */
                f = getR(U, i, l);
                g = -SIGN(sqrt(s), f);
                h = f*g-s;
                setR(U, i, l, f-g);
                
                for ( k = l; k < n; k++ ) {
                    rv1[ k ] = getR(U, i, k)/h;
                } /* end loop k */
                
                for ( j = l; j < m; j++ ) {
                    for ( s = 0.0, k = l; k < n; k++ ) {
                        s += getR(U, j, k) * getR(U, i, k);
                    } /* end loop s */
                    
                    for ( k = l; k < n; k++ ) {
                        setR(U, j, k, getR(U, j, k) + s * rv1[ k ] );
                    } /* end loop k */
                } /* end loop j */
                
                for ( k = l; k < n; k++ ) {
                    setR(U, i, k, getR(U, i, k)*scale);
                } /* end loop k */
            } /* end if ( scale ... */
        } /* end if ( i < m ... */
        anorm = FMAX(anorm, (fabs(wr[i]) + fabs(rv1[i])));
    } /* end loop i */
    
    /* Accumulation of right-hand transformations. */
    for ( i = (int)n-1; i >= 0; i-- ) {
        if ( i < (int)n - 1 ) {
            if ( g != 0 ) {
                for ( j = l; j < (int)n; j++ )  {/* Doulbe divisio to avoid possible underflow. */
                    setR(V, j, i, (getR(U, i, j)/getR(U, i, l))/g);
                } /* end loop j */
                for ( j = l; j < (int)n; j++ ) {
                    for ( s = 0.0, k = l; k < (int)n; k++ ) {
                        s += getR(U, i, k)*getR(V, k, j);
                    } /* end loop s */
                    for ( k = l; k < n; k++ ) {
                        setR(V, k, j, getR(V, k, j) + s*getR(V, k, i));
                    } /* end loop k */
                } /* end loop j */
            } /* end if ( g != 0 ) */
            
            for ( j = l; j < (int)n; j++ ) {
                setR(V, i, j, 0.0);
                setR(V, j, i, 0.0);
            } /* end loop j */
        } /* end if ( i < n - 1 ) */
        
        setR(V, i, i, 1.0);
        g = rv1[i];
        l = i;
    } /* end loop i */
    
    /* Accumulation of left-hand transformations. */
    for ( i = (int)n - 1; i >= 0; i-- ) {
        l = i + 1;
        g = wr[i];
        
        for ( j = l; j < n; j++ ) {
            setR(U, i, j, 0.0);
        } /* end loop j */
        
        if (g) {
            g = 1.0/g;
            
            for ( j = l; j < (int)n; j++ ) {
                for ( s = 0.0, k = l; k < (int)m; k++ ) {
                    s += getR(U, k, i) * getR(U, k, j);
                } /* end loop k */
                f = ( s/getR(U, i, i))*g;
                for ( k = i; k <(int) m; k++ ) {
                    setR(U, k, j, getR(U, k, j) + f * getR(U, k, i));
                } /* end loop k */
            } /* end loop j */
            
            for ( j = i; j < (int)m; j++ ) {
                setR(U, j, i, getR(U, j, i)*g);
            } /* end loop j */
        } /* end if (g) */
        else {
            for ( j = i; j < (int)m; j++ ){
                setR(U, j, i, 0.0);
            } /* end loop j */
        } /* end else */
        
        setR(U, i, i, getR(U, i, i)+1);
    } /* end loop i */
    
    /* Diagonalization of the bidiagonal form : loop over */
    for ( k = (int)n - 1; k >= 0; k-- ) {
        /* Singular values, and over allowed iterations. */
        for ( its = 1; its <= 30; its++ ) {
            flag = 1;
            /* Test for splitting. */
            for ( l = k; l >= 0; l-- ) {
                nm = l - 1;     /* Note that rv1[1] is always zero.*/
                if ((double)(fabs(rv1[l])+anorm) == anorm ) {
                    flag = 0;
                    break;
                } /* end if ((double... */
                if ((double)(fabs(wr[nm]) + anorm) == anorm) {
                    break;
                }
            } /* end loop l */
            
            if (flag ) {
                c = 0.0;
                s = 1.0;
                
                for ( i = l; i <= k; i++ ) {
                    f = s*rv1[i];
                    rv1[i] = c*rv1[i];
                    if ((double)(fabs(f)+anorm) == anorm ) break;
                    g = wr[ i ];
                    h = pythag(f, g);
                    wr[ i ] = h;
                    h = 1.0/h;
                    c = g*h;
                    s = -f*h;
                    for ( j = 1; j <(int) m; j++ ) {
                        y = getR(U, j, nm);
                        z = getR(U, j, i);
                        setR(U, j, nm, y*c+z*s);
                        setR(U, j, i, z*c-y*s);
                    }/* end loop j */
                } /* end loop i */
            } /* end if (flag) */
            
            z = wr[ k ];
            
            if ( l == k ) { /* Convergence */
                if ( z < 0.0 ) {    /* Singular value is made nonegative. */
                    wr[ k ] = -z;
                    for ( j = 0; j < (int)n; j++ ) {
                        setR(V, j, k, -getR(V, j, k));
                    } /* end loop j */
                } /* end if ( z < 0.0 ) */
                break;
            } /* end if (l == k )*/
            
            if ( its == 30 ) {
                mexErrMsgTxt("no convergence in 30 svdcomp iterations");
            }
            x = wr[ l ]; /* Shift from bottom 2-by-2 minor.*/
            nm = k - 1;
            y = wr[nm];
            g = rv1[nm];
            h = rv1[k];
            f = ((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
            g = pythag(f, 1.0);
            f = ((x-z)*(x+z)+h*((y/(f+SIGN(g, f)))-h))/x;
            c=s=1.0;        /* Next QR transformation */
            
            for ( j = l; j <= nm; j++ ) {
                i = j + 1;
                g = rv1[i];
                y = wr[ i ];
                h = s*g;
                g = c*g;
                z = pythag(f, h);
                rv1[j] = z;
                c = f/z;
                s = h/z;
                f = x*c+g*s;
                g = g*c-x*s;
                h = y*s;
                y *= c;
                for ( jj = 0; jj < (int)n; jj++ ) {
                    x = getR(V, jj, j);
                    z = getR(V, jj, i);
                    setR(V, jj, j, x*c+z*s);
                    setR(V, jj, i, z*c-x*s);
                } /* end loop jj */
                z = pythag(f, h);
                
                wr[j] = z;  /* Rotation can be arbitrary if z = 0 */
                if (z) {
                    z = 1.0/z;
                    c = f*z;
                    s = h*z;
                } /* end if (z) */
                f = c*g+s*y;
                x = c*y-s*g;
                
                for ( jj = 0; jj < (int)m; jj++ ) {
                    y = getR(U, jj, j);
                    z = getR(U, jj, i);
                    setR(U, jj, j, y*c+z*s);
                    setR(U, jj, i, z*c-y*s);
                } /* end loop jj */
            } /* end loop j */
            rv1[l] = 0.0;
            rv1[k] = f;
            wr[k] = x;
        } /* end loop its */
    } /* end loop k */
    
    mxFree(rv1);
}

double cmax( const mxArray *v) {
    size_t m,n;
    int i;
    double *vpr;
    double max;
    
    m = ( mxGetM(v) > mxGetN(v) ? mxGetM(v):mxGetN(v) );
    n = ( mxGetM(v) < mxGetN(v) ? mxGetM(v):mxGetN(v) );
    
    if ( (int)n > 1 ) {
        mexErrMsgTxt("max can be used with a vector");
    }
    
    vpr = mxGetPr(v);
    
    max = vpr[0];
    
    if ((int)m == 1 ) {
        return max;
    }
    
    
    for ( i = 1; i < (int)m; i++ ) {
        if ( vpr[i] > max ) max = vpr[i];
    }
    
    return max;
    
}

/* Pseudo Inverse */
mxArray* pinv( const mxArray *A ) {
    size_t m,n;
    double tol, r;
    int i;
    double *wpr;
    mxArray *Pinv;
    mxArray *wsubidiag;
    mxArray *Pinv4;
    mxArray *UTsub;
    mxArray *Pinv1;
    mxArray *UT;
    
    m = mxGetM(A);
    n = mxGetN(A);
    if ( n > m ) {
        mxArray *AT;
        
        mxArray *Pinv2;
        mxArray *Pinv3;
        
        AT = transpose(A);
        Pinv2 = pinv(AT);
        mxDestroyArray(AT);
        
        Pinv3 = transpose(Pinv2);
        mxDestroyArray(Pinv2);
        return Pinv3;
    }
    else {
        mxArray *U = mxCreateDoubleMatrix(m, n, mxREAL);
        mxArray *V = mxCreateDoubleMatrix(n, n, mxREAL);
        mxArray *w = mxCreateDoubleMatrix(n, 1, mxREAL);
        svdcmp(A, U, V, w);
        
        tol = m*(1e-16)*cmax(w);
        
        r = 0;
        
        wpr = mxGetPr(w);
        for ( i = 0; i < n; i++ ) {
            if ( wpr[ i ] > tol ) r += 1;
        }
        
        if ( r < n ) {
            mxArray *Usub;
            mxArray *Vsub;
            mxArray *wsub;
            
            Usub = sub(U, 0, m-1, 0, r-1);
            Vsub = sub(V, 0, n-1, 0, r-1);
            wsub = sub(w, 0, r-1, 0, 0);
            
            
            wsubidiag = idiag(wsub);
            
            Pinv4 = times(Vsub, wsubidiag);
            mxDestroyArray(Vsub);
            mxDestroyArray(wsub);
            mxDestroyArray(wsubidiag);
            
            UTsub = transpose(Usub);
            mxDestroyArray(Usub);
            Pinv = times(Pinv4, UTsub);
            mxDestroyArray(Pinv4);
            mxDestroyArray(UTsub);
            return Pinv;
            
        }
        else {
            mxArray *widiag;
            widiag = idiag(w);
            
            Pinv1 = times(V, widiag);
            mxDestroyArray(V);
            mxDestroyArray(w);
            mxDestroyArray(widiag);
            
            UT = transpose(U);
            mxDestroyArray(U);
            Pinv = times(Pinv1, UT);
            mxDestroyArray(Pinv1);
            mxDestroyArray(UT);
            return Pinv;
        }
    } /* end else */
}


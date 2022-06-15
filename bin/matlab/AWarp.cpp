

#include "mex.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#define inf 1e19;
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))

double UBCases(double a, double b, char c)
{
    double v = 0;
    if( a > 0  && b > 0 )
        v = (a-b)*(a-b);
    else if ( a> 0 && b < 0 )
    {
        if( c=='l' )
            v = a*a;
        
        else
            
            v = (-b)*a*a;
    
    }
    else if ( a < 0 && b > 0 )
    {
        if( c == 't' )
            v = b*b;
        else
            v = (-a)*b*b;
    
    }
    else
        v = 0;

    return v;
}

double dtw_G(double *s, double *t, int c, int ns, int nt )
{
    double d=0;
    double ** D;
    int i,j;
    
    double a1, a2, a3;

    D = (double **) malloc((ns+1)*sizeof(double *));
    for( i = 0 ; i <= ns ; i++ )
        D[i] =(double *)malloc((nt+1)*sizeof(double));
    
    
    
    
    for(i=0;i<ns+1;i++)
    {
       D[i][0] = inf;
      
    }
    
    for(i=0;i<nt+1;i++)
    {
       D[0][i] = inf;
      
    }

    D[0][0]=0;
    
    
    for(i=0;i<ns;i++)
    {
  
        for(j=0;j<nt;j++)
        {
            
            a1 = D[i][j] + (s[i]-t[j])*(s[i]-t[j]);
            
            if( i > 0  && j > 0 )
                a1 = D[i][j] + UBCases(s[i],t[j],'d');
            a2 = D[i+1][j] + UBCases(s[i],t[j],'t');
            a3 = D[i][j+1] + UBCases(s[i],t[j],'l');
            
            D[i+1][j+1]= MIN(a3,MIN(a1,a2));
        }
       
    }
   
    d = sqrt(D[ns][nt]);

    for( i = 0 ; i <= ns ; i++ )
        free(D[i]);
    free(D);

    return d;//(double)matchLength;
}


void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    double *s,*t;
    char  c;
    int ns,nt;
    double *dp;
    
    
    if(nrhs!=2&&nrhs!=3)
    {
        mexErrMsgIdAndTxt( "MATLAB:dtw_c:invalidNumInputs",
                "Two or three inputs required.");
    }
    if(nlhs>1)
    {
        mexErrMsgIdAndTxt( "MATLAB:dtw_c:invalidNumOutputs",
                "dtw_c: One output required.");
    }
    
    
    if(nrhs==2)
    {
        c='u';
    }
    else if(nrhs==3)
    {
        if( !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) ||
                mxGetN(prhs[2])*mxGetM(prhs[2])!=1 )
        {
            mexErrMsgIdAndTxt( "MATLAB:dtw_c:wNotScalar",
                    "dtw: Input W must be a scalar.");
        }
        
       
        c = (int) mxGetScalar(prhs[2]);
    }
    
    
    
    s = mxGetPr(prhs[0]);
    
    t = mxGetPr(prhs[1]);
    
    ns = mxGetM(prhs[0])*mxGetN(prhs[0]);
    
    nt = mxGetM(prhs[1])*mxGetN(prhs[1]);
    
    plhs[0] = mxCreateDoubleMatrix( 1, 2, mxREAL);
    
    dp = mxGetPr(plhs[0]);
    
    dp[0]=dtw_G(s,t,c,ns,nt);
    dp[1] = dp[0];
    
    return;
    
}

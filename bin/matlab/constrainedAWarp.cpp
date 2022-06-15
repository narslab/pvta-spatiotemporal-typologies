
#include "mex.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#define inf 1e19;
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))

double UBCases(double a, double b, char c, int w, int gap)
{
    double v = 0;
    if( a > 0  && b > 0 && gap<= w )
		v = (a-b)*(a-b);
		
	else if (a < 0 && b < 0)
	{
		v = 0;
	}else{
		if (c == 'd'){
			if (a > 0 && b < 0){
				v = (-b)*a*a;
			}else if (a < 0 && b > 0){
				v = (-a)*b*b;
				
			}else{
				v = inf;
			}
		}else if (c == 'l'){
			if (a > 0 && b < 0 && gap <=w){
				v = (-b)*a*a;
			}else if (a < 0 && b > 0){
				v = b*b;
			}else
				v = inf;
			
		}else if (c == 't'){
			if (a > 0 && b < 0 ){
				v = a*a;
			}else if (a < 0 && b > 0 && gap <=w){
				v = (-a)*b*b;
			}else
				v = inf;
		} 
	}
		
	return v;
}


double dtw_G(double *s, double *t, int w, int ns, int nt )
{
	
    double d=0;
    double ** D;
    int i,j;
    double a1, a2, a3;
	double * x;
	double * y;
	int * tx;
	int * ty;
	
	
	
	x = (double *)malloc((ns+1)*sizeof(double));
	y = (double *)malloc((nt+1)*sizeof(double));
	
	tx = (int *)malloc((ns+1)*sizeof(int));
	ty = (int *)malloc((nt+1)*sizeof(int));
	
    D = (double **) malloc((ns+1)*sizeof(double *));
    for( i = 0 ; i <= ns ; i++ )
        D[i] =(double *)malloc((nt+1)*sizeof(double));
    
    if (w < abs(ns-nt)){
		printf("w need to be greater than %d.",abs(ns-nt));
		return(-1);
	}
	
	
	//transfer s to x and append one to them
    for(i=0;i<ns;i++)
    {
       x[i] = s[i];
    }
    x[ns] = 1;
    
    //transfer t to y and append one to them
    for(i=0;i<nt;i++)
    {
       y[i] = t[i];
    }
    y[nt] = 1;
    
    
    //calculate the timestamp of the event in s
    int iit = 0;
    for(i=0;i<ns+1;i++){
       if (x[i] > 0){
		   iit = iit+1;
		   tx[i] = iit;
	   }else{
		   iit = iit + abs(x[i]);
		   tx[i] = iit;		   
	   }
    }
    
    iit = 0;
    for(i=0;i<nt+1;i++){
       if (y[i] > 0){
		   iit = iit+1;
		   ty[i] = iit;
	   }else{
		   iit = iit + abs(y[i]);
		   ty[i] = iit;		   
	   }
    }
	
    for(i=0;i<ns+1;i++)
    {
       D[i][0] = inf;
      
    }
    
    for(i=0;i<nt+1;i++)
    {
       D[0][i] = inf;
      
    }

    D[0][0]=0;
    
    int gap = 0;
    for(i=0;i<ns;i++)
    {
  
        for(j=0;j<nt;j++)
        {
			gap = abs(tx[i] - ty[j]);
            
            if (gap > w && ((j > 0 && ty[j-1]-tx[i] > w) || (i > 0 && tx[i-1]-ty[j] > w))){
				D[i+1][j+1] = inf;
			}
			else{
				a1 = D[i][j] + (s[i]-t[j])*(s[i]-t[j]);
				
				if( i > 0  && j > 0 )
					a1 = D[i][j] + UBCases(s[i],t[j],'d',w,gap);
				a2 = D[i+1][j] + UBCases(s[i],t[j],'l',w,gap);
				a3 = D[i][j+1] + UBCases(s[i],t[j],'t',w,gap);
				
				D[i+1][j+1]= MIN(a3,MIN(a1,a2));
			}
        }
       
    }
      //print the D matrix  
    /*for( i = 0 ; i <= ns ; i++ ){
        for (j = 0 ; j <= nt ; j++){
			
			
			if (D[i][j] >= 10000000000000000000.000000){
					printf("inf ");
			}else{
				printf("%f ", D[i][j] );
			}
		}
        printf("\n");
	}*/
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



/*int main(){
	
    
    
    int c = (int)'u';
    int ns = 8;
    int nt = 5;
    int w = 7;
    int i;
    double s[8] = {1.0, -4.0, 1.0, -4.0, 1.0, 1.0, -4.0, 1.0};
    double t[5] = {1.0, -4.0, 1.0, -4.0, 1.0};
    
    double res = 0;
    
	res = dtw_G(s, t, w, ns, nt);
	printf("%f\n" , res);
	return 0;
	}*/

/*Incomplete LDL factorization using Robust Incomplete Factorization*/

#include "mex.h"
#include "matrix.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>


double signfun(double x)
{
    return (x >= 0) - (x < 0);
}

double absval(double x)
{
    return x * signfun(x);
}


int eyeinit(mwIndex n, mwIndex *zcolstart, mwIndex *zcolend, 
        mwIndex *zrow, double *zval, mwIndex gap)
{
    mwIndex j, next;
    
    for(j = 0;j < n;j++)
    {
        if(j > 0)next = next + j + 1;
        else
            next = 0;
        zcolstart[j] = next;
        zcolend[j] = next;
        zrow[next] = j;
        zval[next] = 1;
    }
/*    for(j = 0;j < gap && j < n;j++)
    {
        if(j > 0)next = next + j + 1;
        else
            next = 0;
        zcolstart[j] = next;
        zcolend[j] = next;
        zrow[next] = j;
        zval[next] = 1;
    }
    
    for(j = gap;j < n;j++)
    {
        next = next + gap;
        zcolstart[j] = next;
        zcolend[j] = next;
        zrow[next] = j;
        zval[next] = 1;
    }
*/
    return 0;
}


/*Compute Matrix Vector Multiplication*/
int matxvec(mwIndex n,mwIndex *ajc,mwIndex *air,
        double *apr,mwIndex *zcolstart,mwIndex *zcolend,
        mwIndex *zrow,double *zval,mwIndex j,double *tmpAzj)
{
    mwIndex col,iter1,iter2,row;
    double val;    
    for(iter1 = zcolstart[j];iter1<=zcolend[j];iter1++)
    {
        col = zrow[iter1];
        val = zval[iter1];
        for(iter2 = ajc[col];iter2 < ajc[col+1];iter2++)
        {
            row = air[iter2];
            tmpAzj[row] += apr[iter2] * val;
        }
    }
    return 0;
}

int addspvecs(mwIndex *zcolstart,mwIndex *zcolend,mwIndex *zrow,
        double *zval,mwIndex i,mwIndex j, double lambda,double eta2)
{
    mwIndex row_coli,row_colj;
    mwIndex iter_coli,iter_colj,tmpind;
    int iter;
    double *sumvecval, zthresh=0;
    if(i<j)
    {
        printf("i must be greater than j\n");
        return 1;
    }    
    sumvecval = (double *)mxCalloc(i+1,sizeof(double));    
    for(iter_coli = zcolstart[i];iter_coli<=zcolend[i];iter_coli++)
    {
        row_coli = zrow[iter_coli];
        sumvecval[row_coli] = zval[iter_coli];
    }
    for(iter_colj=zcolstart[j];iter_colj<=zcolend[j];iter_colj++)
    {
        row_colj = zrow[iter_colj];
        sumvecval[row_colj] += lambda * zval[iter_colj];
    }
    
    for(iter = i;iter >= 0;iter--)zthresh += absval(sumvecval[iter]);
    zthresh = eta2 * zthresh;
    
    tmpind = zcolend[i];
    for(iter=i;iter>=0;iter--)
    {
        if(absval(sumvecval[iter]) > zthresh)
        {
            
            zrow[tmpind] = iter;
            zval[tmpind] = sumvecval[iter];
            tmpind--;
        }
        if(tmpind <= zcolend[i-1])break;
    }
    zcolstart[i] = tmpind + 1;
    mxFree(sumvecval);
    return 0;
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, 
        const mxArray *prhs[])
{   
    /*declaration of variables*/
    /*iteration variables*/
    mwIndex i,j,k,numofnz;
    
    /*Auxillary variables*/
    double tmplij, tmpljj;
    mwIndex tmprowind;

    /*Step 2:*/
    mwSize m, n;/*m:number of rows, n:number of columns*/
    mwSize nzmax;/*nnz of input*/
    
    mwIndex rnzmax;/*rnzmax should be comparable with nzmax*/
    
    mwIndex *ajc, *air;
    double *apr;/*CSC format of prhs[0]*/

    double mu, eta1, eta2;/* mu is the shift, and eta1 is threshold*/
    double *threshold;
    
    /*Step 3: */
    mwIndex *zcolstart, *zcolend, *zrow;
    double *zval;
    mwIndex gap = 5000;
    
    /*Step 4:*/
    double pjj,pij,lambda,zji;
    double  *tmpAzj;
    
    /*Output:*/    
    mwIndex *ljc, *lir;
    double *lpr;
    /*---------------------------------------------*/
    
    /* Step 1: Check input---------*/
    if(nrhs != 4 && nrhs !=5)
    {
        mexErrMsgTxt("Four or Five Input Arguments are required");
    }
    
    if(nlhs != 1)
        mexErrMsgTxt("One Output Argument Required\n");
        
    /*Check whether prhs[] is sparse*/
    if(mxIsSparse(prhs[0]) != 1)
        mexErrMsgTxt("Input Matrix Must Be Sparse\n");
    /* ------------------------------*/
    
    
    /* Step 2: Get Input Values------------*/
    /*Read m and n from prhs[0]*/
    m = mxGetM(prhs[0]);
    n = mxGetN(prhs[0]);
    nzmax = mxGetNzmax(prhs[0]);
    
    /*CSC format of A=prhs[0]*/
    ajc = mxGetJc(prhs[0]);
    air = mxGetIr(prhs[0]);
    apr = mxGetPr(prhs[0]);
       
    
    /* Get shift*/
    mu = mxGetScalar(prhs[1]);
    
    /*Get threshold*/
    eta1 = mxGetScalar(prhs[2]);
    threshold = (double *)mxCalloc(n,sizeof(double)); /* why are you doing it?*/
    
    for(j=0;j<n;j++)
    {
        for(i=ajc[j];i<ajc[j+1];i++)threshold[j] += absval(apr[i]);
        threshold[j] = eta1 * threshold[j];
    }
    
    /*Get threshold parameter for Z*/
    eta2 = mxGetScalar(prhs[3]);
    
    /* Get allocation parameter*/
    if(nrhs == 5)
    {
        gap = (mwIndex)mxGetScalar(prhs[4]);
        if(gap > n)gap = n;
    }
    
    rnzmax = n*n; /*lingfei*/
    if(gap * m > rnzmax) rnzmax = (mwIndex)(gap*m); /*why are you doing it? */
            
    
    /*---------------------------------------*/
    

    /*Step 3: Initialization of Z and L----------    */
    zcolstart = (mwIndex *)mxCalloc(n, sizeof(mwIndex));
    zcolend = (mwIndex *)mxCalloc(n, sizeof(mwIndex));
    zrow = (mwIndex *)mxCalloc(rnzmax,sizeof(mwIndex));
    zval = (double *)mxCalloc(rnzmax,sizeof(double));
    
    if(zrow == NULL || zval == NULL){
        printf("Out of memory, please use a smaller gap value\n");
        return;
    }
    
    eyeinit(n, zcolstart, zcolend, zrow, zval, gap);/* Let Z be eye(n,n)*/
    
    /*---------------------------------------*/
    
    
    /* Step : Output    */
    plhs[0] = mxCreateSparse(n,n,rnzmax,mxREAL);
    
    ljc = mxGetJc(plhs[0]);
    lir = mxGetIr(plhs[0]);
    lpr = mxGetPr(plhs[0]);
    
    
    
    /*Step 4: Compute L*/
    
    numofnz = 0;
    
    for(j=0;j<n;j++)
    {
        
        pjj = 0;
        tmpAzj = (double *)mxCalloc(m,sizeof(double));
        matxvec(n,ajc,air,apr,zcolstart,zcolend,zrow,zval,j,tmpAzj);
        for(k=0;k<m;k++)
            if(tmpAzj[k] != 0)pjj += tmpAzj[k] * tmpAzj[k];
        if(mu != 0){
            for(i = zcolstart[j]; i <= zcolend[j]; i++)
                pjj -= mu * mu * zval[i] * zval[i];
        }
        ljc[j] = numofnz;
        lir[numofnz] = j;
        lpr[numofnz] = sqrt(absval(pjj));
        tmpljj = lpr[numofnz];
        numofnz = numofnz + 1;
        if(tmpljj < 2.2e-16 ){
            lpr[numofnz-1] = (threshold[j]>0)?threshold[j]:2.2e-16;
            continue;
        }
        for(i=j+1;i<n;i++)
        {
            pij = 0;
            for(k=ajc[i];k<ajc[i+1];k++)
            {
                tmprowind = air[k];
            	pij += apr[k] * tmpAzj[tmprowind];
            }
            zji = 0;
            for(k=zcolend[i];k>=zcolstart[i];k--)
            {
                if(zrow[k]==j)zji = zval[k];
            }
            pij -= mu * mu * zji;

            lambda = pij / pjj;
            tmplij = pij / tmpljj * signfun(pjj);
            if(absval(tmplij) > threshold[i])
            {
                lir[numofnz] = i;
                lpr[numofnz] = tmplij;
                numofnz = numofnz + 1;
                addspvecs(zcolstart,zcolend,zrow,zval,i,j,-lambda,eta2);
            }
        }
        mxFree(tmpAzj);
    }
    ljc[n] = numofnz;
}

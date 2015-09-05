/*******************************************************************************
 *   PRIMME PReconditioned Iterative MultiMethod Eigensolver
 *   Copyright (C) 2015 College of William & Mary,
 *   James R. McCombs, Eloy Romero Alcalde, Andreas Stathopoulos, Lingfei Wu
 *
 *   This file is part of PRIMME.
 *
 *   PRIMME is free software; you can redistribute it and/or
 *   modify it under the terms of the GNU Lesser General Public
 *   License as published by the Free Software Foundation; either
 *   version 2.1 of the License, or (at your option) any later version.
 *
 *   PRIMME is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *   Lesser General Public License for more details.
 *
 *   You should have received a copy of the GNU Lesser General Public
 *   License along with this library; if not, write to the Free Software
 *   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 *******************************************************************************
 * File: PRIMME_mex.c
 * 
 * Purpose - PRIMME MEX interface.
 * 
 * Currently, we only support sequential computing. The PRIMME MEX
 * reads the inputs from MATLAB, constructs the necessary structures
 * and then calls PRIMME. The desired results are then returned to MATLAB.
 * 
 * Matrix-vector and preconditioning functions are performed by callbacks
 * to MATLAB functions.
 *
 * We provide different levels of function calls (similar to MATLAB eigs())
 *
 * Output: [evals, evecs, norms, primmeout]
 * 
 * Function call:  PRIMME_eigs(flag, dim)
 *                 PRIMME_eigs(flag, dim, numEvals)
 *                 PRIMME_eigs(flag, dim, numEvals, target)
 *                 PRIMME_eigs(flag, dim, numEvals, target, method)
 *                 PRIMME_eigs(flag, dim, numEvals, target, method, opts)
 *
 * Input: [dim, numEvals, target, method, opts]
 *
 *[1] 'flag' is to mark the input matrix as real or complex
 *[2] 'dim' is the dimension of the symmetric/heritian matrix
 *[3] 'numEvals' is number of eigenvalues returned               [default = 1]
 *[4] 'target' is largest, smallest, closest_geq, closest_leq, and 
 *     closest_abs                                 [default = primme_smallest]
 *[5] 'Method' is DYNAMIC, DEFAULT_MIN_TIME, DEFAULT_MIN_MATVECS, Arnoldi,
 *     GD, .. and so on                                    [default = DYNAMIC]
 *[6] 'opts' is an option structure which contain following parameters 
 *   in the primme_params structure: 
 *   (0) opts.aNorm: the estimate norm value of matrix A         [default = 0]
 *   (1) opts.eps: desired computing accuracy                [default = 1e-12]
 *   (2) opts.numTargetShifts: shifts for interior eigenvalues   [default = 0]
 *   (3) opts.targetShifts: pointer to get each shift for interior eigenvalues
 *   (4) opts.initSize: initial guesses/constraints              [default = 0]
 *   (5) opts.numOrthoConst:                                     [default = 0]
 *   (6) opts.locking: 0 or 1 
 *   (7) opts.dynamicMethodSwitch: from -3 to 1
 *   (8) opts.maxBasisSize: maximum basis size allowed in the main iteration
 *   (9) opts.minRestartSize: minimum Ritz vectors to restart
 *  (10) opts.maxBlockSize:                                      [default = 1]
 *  (11) opts.maxMatvecs:                                  [default = INT_MAX]
 *  (12) opts.maxOuterIterations:                          [default = INT_MAX]
 *  (13) opts.restartingParams.scheme:                [default = primme_thick]
 *  (14) opts.restartingParams.maxPrevRetain:                    [default = 1]
 *  (15) opts.precondition: set to 1 if preconditioning is to be performed,
 *       make sure the applyPreconditioner is not NULL           [default = 0]
 *  (16) opts.robustShifts: set to 1 to use robustShifting
 *  (17) opts.maxInnerIterations: 0 (no inner iterations), = k (perform 
 *       at most k inner iterations per outer step)
 *  (18) opts.LeftQ: 0 or 1
 *  (19) opts.LeftX: 0 or 1
 *  (20) opts.RightQ: 0 or 1
 *  (21) opts.RightX: 0 or 1
 *  (22) opts.SkewQ: 0 or 1
 *  (23) opts.SkewX: 0 or 1
 *  (24) opts.relTolBase: a legacy from classical JDQR
 *  (25) opts.convTest: how to stop the inner QMR method
 *  (26) opts.printLevel: 0-5 (different level reporting)         [default =0]
 *  (27) opts.outputFileName: support two output file name
 *  (28) opts.iseed: set iseed value for initialization
 *  (29) opts.intWorkSize: memory size for int workspace
 *  (30) opts.realWorkSize: memory size for real workspace  
 *
 *  For details about PRIMME parameters, methods, and settings see ../readme.txt
 *
 ******************************************************************************/

#include "mex.h"
#include "lapack.h"
#include "blas.h"
#include "primme.h" 
#include "wtime.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#ifndef max
#define max(a, b) ((a) > (b) ? (a) : (b))
#endif

void MatrixMatvec_d(void *x, void *y, int *blockSize, primme_params *primme);
void MatrixMatvec_z(void *x, void *y, int *blockSize, primme_params *primme);
void Preconditioner_d(void *x, void *y, int *blockSize, primme_params *primme);
void Preconditioner_z(void *x, void *y, int *blockSize, primme_params *primme);

double Matvec_mex_timer = 0.0L;
char *outputfilename;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   double *evals, *rnorms;
   void *EVecs;
   primme_params primme;
   mwSize i, j;
   mwSize ndim;
   mxArray *tmp; /* tmp stores each field value pointer of opts*/
   char *read_target_string = NULL;
   char *read_method_string = NULL; /* will be used in future version*/   
   char *read_projection_string = NULL; 
   double *read_iseed = NULL;
   double *read_initialevecs = NULL;
   double *read_initialevecsimag = NULL;

   /* Timing vars */
   double wt1,wt2;

   /* check: The number of input arguments are between 1 and 6 */
   if (nrhs == 0)
      mexErrMsgTxt("Must have at least one input arguments");
   else if (nrhs > 6)
      mexErrMsgTxt("Must have at most 6 input arguments");

   /* check: The number of output arguments are between 1 and 4 */
   if (nlhs == 0)
      mexErrMsgTxt("Must have at least one output arguments");
   else if (nlhs > 4)
      mexErrMsgTxt("Must have at most 4 output arguments");

   mexPrintf("input and output arguments check finished\n");

   primme_initialize(&primme);
   primme_preset_method method = DYNAMIC;

   mexPrintf("initialization finished...\n");

   /* check: set different initial parameters according to the number of input arguments*/
   if (nrhs >= 1)
   {
      if (mxIsEmpty(prhs[0])) {
         mexErrMsgTxt("Input matrix cannot be empty matrix");
      }
   }
   if (nrhs >= 2) {
      if (!mxIsEmpty(prhs[1])) {
         primme.n = (int)mxGetScalar(prhs[1]); /*get the dimension of the matrix*/
         primme.nLocal = primme.n;  /*set for sequential programs */
      }
      else
         mexErrMsgTxt("The dimension of the constructed matrix is not provided");
      mexPrintf("primme.n is %d\n", primme.n);
   }
   if (nrhs >= 3) {
      if ((int)mxGetScalar(prhs[2])<= primme.n)
         primme.numEvals = (int)mxGetScalar(prhs[2]); /*get the number of eigenvalues*/
      else
         mexErrMsgTxt("The number of eigenvalues must be less than the dimension of the matrix");
      mexPrintf("Number of eigenvalues is %d\n", primme.numEvals);
   }
   if (nrhs >= 4) {
      if (!mxIsChar(prhs[3]) || (mxGetM(prhs[3]) != 1 ))
         mexErrMsgTxt("The target argument must be a string");
      /* primme_smallest    --- 0 --- "SA"
       * primme_largest     --- 1 --- "LA"
       * primme_closest_geq --- 2 --- "CGT"
       * primme_closest_leq --- 3 --- "CLT"
       * primme_closest_abs --- 4 --- "CT" */
      else {
         read_target_string = mxArrayToString(prhs[3]);
         if (!strcmp(read_target_string,"SA")) {
            primme.target = 0; /*get the target method to find eigenvalues*/
            mexPrintf("primme.target is primme_smallest\n");
         }
         if (!strcmp(read_target_string,"LA")) {
            primme.target = 1;
            mexPrintf("primme.target is primme_largest\n");
         }
         if (!strcmp(read_target_string,"CGT")) {
            primme.target = 2;
            mexPrintf("primme.target is primme_closest_geq\n");
         }
         if (!strcmp(read_target_string,"CLT")) {
            primme.target = 3;
            mexPrintf("primme.target is primme_closest_leq\n");
         }
         if (!strcmp(read_target_string,"CT")) {
            primme.target = 4;
            mexPrintf("primme.target is primme_closest_abs\n");
         }
      }
      mxFree(read_target_string);
   }
   if (nrhs >= 5) {
      if ((int)mxGetScalar(prhs[4]) >=0 && (int)mxGetScalar(prhs[4])<=14)
         method = (int)mxGetScalar(prhs[4]); /*get the solver method the user chooses*/ 
      else
         mexErrMsgTxt("The method number must be in 0-14");  
      /*method   One of the following 12 enum methods:
       * typedef enum{
       * DYNAMIC,                  ---0: Switches dynamically to the best method
       * DEFAULT_MIN_TIME,         ---1: Currently set at JDQMR_ETol
       * DEFAULT_MIN_MATVECS,      ---2: Currently set at GD+block
       * Arnoldi,                  ---3: obviously not an efficient choice 
       * GD,                       ---4: classical block Generalized Davidson 
       * GD_plusK,                 ---5: GD+k block GD with recurrence restarting
       * GD_Olsen_plusK,           ---6: GD+k with approximate Olsen precond.
       * JD_Olsen_plusK,           ---7: GD+k, exact Olsen (two precond per step)
       * RQI,                      ---8: Rayleigh Quotient Iteration. Also INVIT,
       *                                 but for INVIT provide targetShifts
       * JDQR,                     ---9: Original block, Jacobi Davidson
       * JDQMR,                    ---10: Our block JDQMR method (similar to JDCG)
       * JDQMR_ETol,               ---11: Slight, but efficient JDQMR modification
       * SUBSPACE_ITERATION,       ---12: equiv. to GD(block,2*block)
       * LOBPCG_OrthoBasis,        ---13: equiv. to GD(nev,3*nev)+nev
       * LOBPCG_OrthoBasis_Window  ---14: equiv. to GD(block,3*block)+block nev>block 
       * }*/
      switch(method) {
         case 0:
            mexPrintf("The solver method is DYNAMIC\n");
            break;
         case 1:
            mexPrintf("The solver method is DEFAULT_MIN_TIME\n");
            break;
         case 2:
            mexPrintf("The solver method is DEFAULT_MIN_MATVECS\n");
            break;
         case 3:                
            mexPrintf("The solver method is Arnoldi\n");
            break;
         case 4:
            mexPrintf("The solver method is GD\n");     
            break;
         case 5:
            mexPrintf("The solver method is GD_plusK\n");   
            break;
         case 6:  
            mexPrintf("GD_Olsen_plusK\n");
            break;
         case 7:     
            mexPrintf("JD_Olsen_plusK\n");   
            break;
         case 8:  
            mexPrintf("The solver method is RQI\n");     
            break;
         case 9:
            mexPrintf("The solver method is JDQR\n");   
            break;
         case 10:  
            mexPrintf("The solver method is JDQMR\n");   
            break;
         case 11:  
            mexPrintf("The solver method is JDQMR_ETol\n");    
            break;
         case 12: 
            mexPrintf("The solver method is SUBSPACE_ITERATION\n");     
            break;
         case 13:
            mexPrintf("The solver method is LOBPCG_OrthoBasis\n"); 
            break;
         case 14:    
            mexPrintf("The solver method is LOBPCG_OrthoBasis_Window\n");   
            break;
      }  
   }
   primme_set_method(method, &primme);
   if (nrhs == 6) {
      if (!mxIsStruct(prhs[5]))
         mexErrMsgTxt("The opts must be a structure");
      else {
         tmp = mxGetField(prhs[5], 0, "aNorm");
         if (tmp != NULL) {
            primme.aNorm = (double)mxGetScalar(tmp);
            mexPrintf("primme.aNorm is:%e\n",primme.aNorm);
         }
         tmp = mxGetField(prhs[5], 0, "eps");
         if (tmp != NULL) {
            primme.eps = (double)mxGetScalar(tmp);
            mexPrintf("primme.eps is:%e\n",primme.eps);
         }
         tmp = mxGetField(prhs[5], 0, "numTargetShifts");
         if (tmp != NULL) {
            primme.numTargetShifts = (int)mxGetScalar(tmp);
            mexPrintf("primme.numTargetShifts is:%d\n",primme.numTargetShifts);
         }
         tmp = mxGetField(prhs[5], 0, "targetShifts");
         if (tmp != NULL) {
            primme.targetShifts = mxGetPr(tmp);   
            for (i=0; i< primme.numTargetShifts; i++)               
               mexPrintf("primme.targetShifts[%d] is:%e\n", i, primme.targetShifts[i]);     
         }
         tmp = mxGetField(prhs[5], 0, "initSize");
         if (tmp != NULL) {
            primme.initSize = (int)mxGetScalar(tmp);
            mexPrintf("primme.initSize is:%d\n",primme.initSize);
         }
         tmp = mxGetField(prhs[5], 0, "numOrthoConst");
         if (tmp != NULL) {
            primme.numOrthoConst = (int)mxGetScalar(tmp);
            mexPrintf("primme.numOrthoConst is:%d\n",primme.numOrthoConst);
         }
         tmp = mxGetField(prhs[5], 0, "locking");
         if (tmp != NULL) {
            primme.locking = (int)mxGetScalar(tmp);  
            mexPrintf("primme.locking is:%d\n",primme.locking);   
         }
         tmp = mxGetField(prhs[5], 0, "dynamicMethodSwitch");
         if (tmp != NULL) {
            primme.dynamicMethodSwitch = (int)mxGetScalar(tmp);  
            mexPrintf("primme.dynamicMethodSwitch is:%d\n",primme.dynamicMethodSwitch);     
         }
         tmp = mxGetField(prhs[5], 0, "maxBasisSize");
         if (tmp != NULL) {
            primme.maxBasisSize = (int)mxGetScalar(tmp);
            mexPrintf("primme.maxBasisSize is:%d\n",primme.maxBasisSize);        
         }
         tmp = mxGetField(prhs[5], 0, "minRestartSize");
         if (tmp != NULL) {
            primme.minRestartSize = (int)mxGetScalar(tmp);   
            mexPrintf("primme.minRestartSize is:%d\n",primme.minRestartSize);      
         }
         tmp = mxGetField(prhs[5], 0, "maxBlockSize");
         if (tmp != NULL) {
            primme.maxBlockSize = (int)mxGetScalar(tmp);  
            mexPrintf("primme.maxBlockSize is:%d\n",primme.maxBlockSize);   
         }
         tmp = mxGetField(prhs[5], 0, "maxMatvecs");
         if (tmp != NULL) {
            primme.maxMatvecs = (int)mxGetScalar(tmp);   
            mexPrintf("primme.maxMatvecs is:%d\n",primme.maxMatvecs);   
         }
         tmp = mxGetField(prhs[5], 0, "maxOuterIterations");
         if (tmp != NULL) {
            primme.maxOuterIterations = (int)mxGetScalar(tmp);   
            mexPrintf("primme.maxOuterIterations is:%d\n",primme.maxOuterIterations);    
         }
         tmp = mxGetField(prhs[5], 0, "scheme");
         if (tmp != NULL) {
            primme.restartingParams.scheme = (int)mxGetScalar(tmp);   
            mexPrintf("primme.restartingParams.scheme is:%d\n",primme.restartingParams.scheme);    
         }
         tmp = mxGetField(prhs[5], 0, "maxPrevRetain");
         if (tmp != NULL) {
            primme.restartingParams.maxPrevRetain = (int)mxGetScalar(tmp);   
            mexPrintf("primme.restartingParams.maxPrevRetain is:%d\n",primme.restartingParams.maxPrevRetain);    
         }
         tmp = mxGetField(prhs[5], 0, "precondition");
         if (tmp != NULL) {
            primme.correctionParams.precondition = (int)mxGetScalar(tmp); 
            mexPrintf("primme.correctionParams.precondition is:%d\n",primme.correctionParams.precondition);      
         }
         tmp = mxGetField(prhs[5], 0, "robustShifts");
         if (tmp != NULL) {
            primme.correctionParams.robustShifts = (int)mxGetScalar(tmp);  
            mexPrintf("primme.correctionParams.robustShifts is:%d\n",primme.correctionParams.robustShifts);   
         }
         tmp = mxGetField(prhs[5], 0, "maxInnerIterations");
         if (tmp != NULL) {
            primme.correctionParams.maxInnerIterations = (int)mxGetScalar(tmp); 
            mexPrintf("primme.correctionParams.maxInnerIterations is:%d\n",primme.correctionParams.maxInnerIterations);   
         }
         tmp = mxGetField(prhs[5], 0, "LeftQ");
         if (tmp != NULL) {
            primme.correctionParams.projectors.LeftQ = (int)mxGetScalar(tmp); 
            mexPrintf("primme.correctionParams.projectors.LeftQ is:%d\n",primme.correctionParams.projectors.LeftQ);    
         }
         tmp = mxGetField(prhs[5], 0, "LeftX");
         if (tmp != NULL) {
            primme.correctionParams.projectors.LeftX = (int)mxGetScalar(tmp);   
            mexPrintf("primme.correctionParams.projectors.LeftX is:%d\n",primme.correctionParams.projectors.LeftX);    
         }
         tmp = mxGetField(prhs[5], 0, "RightQ");
         if (tmp != NULL) {
            primme.correctionParams.projectors.RightQ = (int)mxGetScalar(tmp);  
            mexPrintf("primme.correctionParams.projectors.RightQ is:%d\n",primme.correctionParams.projectors.RightQ);      
         }
         tmp = mxGetField(prhs[5], 0, "RightX");
         if (tmp != NULL) {
            primme.correctionParams.projectors.RightX = (int)mxGetScalar(tmp);  
            mexPrintf("primme.correctionParams.projectors.RightX is:%d\n",primme.correctionParams.projectors.RightX);       
         }
         tmp = mxGetField(prhs[5], 0, "SkewQ");
         if (tmp != NULL) {
            primme.correctionParams.projectors.SkewQ = (int)mxGetScalar(tmp);  
            mexPrintf("primme.correctionParams.projectors.SkewQ is:%d\n",primme.correctionParams.projectors.SkewQ);      
         }
         tmp = mxGetField(prhs[5], 0, "SkewX");
         if (tmp != NULL) {
            primme.correctionParams.projectors.SkewX = (int)mxGetScalar(tmp);   
            mexPrintf("primme.correctionParams.projectors.SkewX is:%d\n",primme.correctionParams.projectors.SkewX);       
         }
         tmp = mxGetField(prhs[5], 0, "relTolBase");
         if (tmp != NULL) {
            primme.correctionParams.relTolBase = (double)mxGetScalar(tmp);   
            mexPrintf("primme.correctionParams.relTolBase is:%e\n",primme.correctionParams.relTolBase);    
         }
         tmp = mxGetField(prhs[5], 0, "convTest");
         if (tmp != NULL) {
            primme.correctionParams.convTest = (int)mxGetScalar(tmp);  
            mexPrintf("primme.correctionParams.convTest is:%d\n",primme.correctionParams.convTest);     
         }
         tmp = mxGetField(prhs[5], 0, "printLevel");
         if (tmp != NULL) {
            primme.printLevel = (int)mxGetScalar(tmp);
            mexPrintf("primme.printLevel is:%d\n",primme.printLevel);     
         }
         tmp = mxGetField(prhs[5], 0, "outputFileName");
         if (tmp != NULL) {  
            outputfilename = mxArrayToString(tmp);
            primme.outputFile = fopen(outputfilename, "a+");
            mexPrintf("primme.outputFileName is:%s\n",outputfilename);
            ndim = mxGetN(prhs[5]); /*get the column dimension of struct*/
         }
         tmp = mxGetField(prhs[5], 0, "iseed");
         if (tmp != NULL) {
            read_iseed = mxGetPr(tmp);   
            for (i=0; i< 5; i++)
               primme.iseed[i] = (int)read_iseed[i];
            mexPrintf("primme.iseed is:[%d, %d, %d, %d]\n",
                  primme.iseed[0], primme.iseed[1], primme.iseed[2], primme.iseed[3]);  
         }
         tmp = mxGetField(prhs[5], 0, "intWorkSize");
         if (tmp != NULL) {
            primme.intWorkSize = (int)mxGetScalar(tmp);   
            mexPrintf("primme.intWorkSize is:%d\n",primme.intWorkSize);     
         }
         tmp = mxGetField(prhs[5], 0, "realWorkSize");
         if (tmp != NULL) {
            primme.realWorkSize = (int)mxGetScalar(tmp);
            mexPrintf("primme.realWorkSize is:%d\n",primme.realWorkSize);                    
         }
         tmp = mxGetField(prhs[5], 0, "initialevecs");
         if (tmp != NULL) {
            read_initialevecs = mxGetPr(tmp);
            read_initialevecsimag = mxGetPi(tmp); 
         }
      }
   }

   mexPrintf("finish reading input parameters...\n");
   mexPrintf("ready to judge if the matrix A is real or complex\n");

   if (!mxIsComplex(prhs[0]))
   {
      mexPrintf("The matrix A is real\n");
      /* set Matrix-vector multiplication function */
      primme.matrixMatvec = MatrixMatvec_d; 
      /* set Preconditioning function */
      primme.applyPreconditioner = Preconditioner_d; 

      double *evecs = (double *)EVecs;
      /* Compute the space needed to allocate for evecs */
      int num;
      /*When locking is 0, num gets the larger of the two parameters*/
      num = max(primme.numEvals,primme.initSize); 
      num = num + primme.numOrthoConst; 
      evals = (double *)mxCalloc(primme.numEvals, sizeof(double));
      evecs = (double *)mxCalloc(primme.n*(num+primme.maxBlockSize), sizeof(double));
      rnorms = (double *)mxCalloc(primme.numEvals, sizeof(double));

      mexPrintf("evals, evecs, rnorms memory allocation succeed\n");

      double *Outevals, *Outevecs, *Outrnorms, *Outstates;

      if (nlhs >= 1) {
         plhs[0] = mxCreateDoubleMatrix(primme.numEvals,1,mxREAL);
         Outevals = mxGetPr(plhs[0]);
      }
      if (nlhs >= 2) {  
         plhs[1] = mxCreateDoubleMatrix(primme.n, num, mxREAL);
         Outevecs = mxGetPr(plhs[1]);
      }
      if (nlhs >= 3) {
         plhs[2] = mxCreateDoubleMatrix(primme.numEvals,1,mxREAL);
         Outrnorms = mxGetPr(plhs[2]);
      }
      if (nlhs >= 4) {
         plhs[3] = mxCreateDoubleMatrix(5,1,mxREAL);
         Outstates = mxGetPr(plhs[3]);
      }     

      /* Show initial configurations */
      if (outputfilename != NULL)
         primme_display_params(primme);

      mexPrintf("\n Starting dprimme .... \n");

      int ret = dprimme(NULL,NULL,NULL,&primme);
      mexPrintf("ret value for reporting memory is %d\n", ret);
      mexPrintf("The intworksize is %d\n", primme.intWorkSize);
      mexPrintf("The realworksize is %d\n", primme.realWorkSize);

      /* ------------------------ */
      /* Initial guess (optional) */
      /* ------------------------ */
      if (read_initialevecs != NULL) {
         for (i =0;i<(primme.initSize+primme.numOrthoConst)*primme.n;i++) {
            evecs[i] = read_initialevecs[i];
         }
      }

      /* ------------- */
      /*  Call primme  */
      /* ------------- */

      wt1 = primme_get_wtime(); 

      ret = dprimme(evals, evecs, rnorms, &primme);

      wt2 = primme_get_wtime();


      if (ret == 0)
         mexPrintf("dprimme return value is %d, success\n", ret);
      else
         mexPrintf("dprimme return value is %d, fail\n", ret);

      mexPrintf("Wallclock Runtime   : %f seconds\n", wt2-wt1);
      mexPrintf("Matvec_MEX Time     : %f seconds\n", Matvec_mex_timer);


      /* Show recommended method for future runs */
      if (outputfilename != NULL) {
         for (i=0; i < primme.numEvals; i++) {
            fprintf(primme.outputFile, "Eval[%lu]: %-22.15E rnorm: %-22.15E\n", i+1,evals[i], rnorms[i]); 
         }
         fprintf(primme.outputFile, " %d eigenpairs converged\n", primme.initSize);
         fprintf(primme.outputFile, "Tolerance : %-22.15E\n", primme.aNorm*primme.eps);
         fprintf(primme.outputFile, "Iterations: %-d\n", primme.stats.numOuterIterations); 
         fprintf(primme.outputFile, "Restarts  : %-d\n", primme.stats.numRestarts);
         fprintf(primme.outputFile, "Matvecs   : %-d\n", primme.stats.numMatvecs);
         fprintf(primme.outputFile, "Preconds  : %-d\n", primme.stats.numPreconds); 
         fprintf(primme.outputFile, "\n#,%d,%.1f\n\n", primme.stats.numMatvecs, wt2-wt1);

         switch (primme.dynamicMethodSwitch) {
            case -1: fprintf(primme.outputFile, "Recommended method for next run: DEFAULT_MIN_MATVECS\n"); break;
            case -2: fprintf(primme.outputFile, "Recommended method for next run: DEFAULT_MIN_TIME\n"); break;
            case -3: fprintf(primme.outputFile, "Recommended method for next run: DYNAMIC (close call)\n"); break;
         }

         fprintf(primme.outputFile, "Wallclock Runtime : %f seconds\n", wt2-wt1);
         fprintf(primme.outputFile, "Matvec_MEX Time   : %f seconds\n", Matvec_mex_timer);
         /*           primme_display_params(primme);*/
         fclose(primme.outputFile);           
      }

      if (nlhs >= 1) {  
         for (i=0;i<primme.numEvals;i++) {
            Outevals[i] = evals[i];
         }
      }

      if (nlhs >= 2) {
         for (i=0; i<(primme.n*num);i++) {
            Outevecs[i] = evecs[i];
         }
      }

      if (nlhs >= 3) {  
         for (i=0;i<primme.numEvals;i++) {
            Outrnorms[i] = rnorms[i];
         }
      }

      if (nlhs >= 4) {
         Outstates[0] = (double)primme.stats.numOuterIterations;
         Outstates[1] = (double)primme.stats.numRestarts;
         Outstates[2] = (double)primme.stats.numMatvecs;
         Outstates[3] = (double)primme.stats.numPreconds;
         Outstates[4] = primme.initSize;
      } 

      mxFree(evals);
      mxFree(evecs);
      mxFree(rnorms);
      primme_Free(&primme);
      if (outputfilename != NULL)
         mxFree(outputfilename);
   }
   else {
      mexPrintf("The matrix A is complex\n");
      /* set the Matrix-vector multiplication */
      primme.matrixMatvec = MatrixMatvec_z; 
      /* set the preconditioning function */
      primme.applyPreconditioner = Preconditioner_z; 
      /*set the method for solving the problem. */   
      primme_set_method(method, &primme); 

      Complex_Z *evecs =(Complex_Z *)EVecs;
      /* Compute the space needed to allocate for evecs */
      int num;
      /*When locking is 0, num gets the larger of the two parameters*/
      num = max(primme.numEvals,primme.initSize); 
      num = num + primme.numOrthoConst; 
      evals = (double *)mxCalloc(primme.numEvals, sizeof(double));
      evecs = (Complex_Z *)mxCalloc(primme.n*(num+primme.maxBlockSize), sizeof(Complex_Z));
      rnorms = (double *)mxCalloc(primme.numEvals, sizeof(double));

      mexPrintf("evals, evecs, rnorms memory allocation succeed\n");

      double *Outevals, *OutevecsR, *OutevecsI, *Outrnorms, *Outstates;

      if (nlhs >= 1) {
         plhs[0] = mxCreateDoubleMatrix(primme.numEvals,1,mxREAL);
         Outevals = mxGetPr(plhs[0]);
      }
      if (nlhs >= 2) {  
         plhs[1] = mxCreateDoubleMatrix(primme.n, num, mxCOMPLEX);
         OutevecsR = mxGetPr(plhs[1]);
         OutevecsI = mxGetPi(plhs[1]);
      }
      if (nlhs >= 3) {
         plhs[2] = mxCreateDoubleMatrix(primme.numEvals,1,mxREAL);
         Outrnorms = mxGetPr(plhs[2]);
      }
      if (nlhs >= 4) {
         plhs[3] = mxCreateDoubleMatrix(4,1,mxREAL);
         Outstates = mxGetPr(plhs[3]);
      }

      /* Show initial configuration */
      if (outputfilename != NULL)
         primme_display_params(primme);

      mexPrintf("\n Starting zprimme .... \n");

      /* compute the needed primme.realWorkSize and primme.intWorkSize */
      int ret = zprimme(NULL,NULL,NULL,&primme); 
      mexPrintf("ret value for reporting memory is %d\n", ret);
      mexPrintf("The intworksize is %d\n", primme.intWorkSize);
      mexPrintf("The realworksize is %d\n", primme.realWorkSize);

      /* ------------------------ */
      /* Initial guess (optional) */
      /* ------------------------ */

      if (read_initialevecs != NULL)
         for (i =0;i<(primme.initSize+primme.numOrthoConst)*primme.n;i++) {
            evecs[i].r = read_initialevecs[i];
            evecs[i].i = read_initialevecsimag[i];
         }

      /* ------------- */
      /*  Call primme  */
      /* ------------- */

      wt1 = primme_get_wtime(); 

      ret = zprimme(evals, evecs, rnorms, &primme);

      wt2 = primme_get_wtime();


      if (ret == 0)
         mexPrintf("zprimme return value is %d, success\n", ret);
      else
         mexPrintf("zprimme return value is %d, fail\n", ret);

      mexPrintf("Wallclock Runtime   : %f seconds\n", wt2-wt1);
      mexPrintf("Matvec_MEX Time     : %f seconds\n", Matvec_mex_timer);


      /* Show recommended method for future runs */
      if (outputfilename != NULL) {  
         switch (primme.dynamicMethodSwitch) {
            case -1: fprintf(primme.outputFile, "Recommended method for next run: DEFAULT_MIN_MATVECS\n"); break;
            case -2: fprintf(primme.outputFile, "Recommended method for next run: DEFAULT_MIN_TIME\n"); break;
            case -3: fprintf(primme.outputFile, "Recommended method for next run: DYNAMIC (close call)\n"); break;
         }

         fprintf(primme.outputFile, "Wallclock Runtime : %f seconds\n", wt2-wt1);
         fprintf(primme.outputFile, "Matvec_MEX Time   : %f seconds\n", Matvec_mex_timer);
         fclose(primme.outputFile);
      }


      if (nlhs >= 1) {
         for (i=0; i<primme.numEvals; i++) {
            Outevals[i] = evals[i];
         }
      }

      if (nlhs >= 2) {
         for (i=0; i<(primme.n*num);i++) {
            OutevecsR[i] = evecs[i].r;
            OutevecsI[i] = evecs[i].i;
         }
      }

      if (nlhs >= 3) {
         for (i=0; i<primme.numEvals; i++) {
            Outrnorms[i] = rnorms[i];
         }
      }

      if (nlhs >= 4) {
         Outstates[0] = (double)primme.stats.numOuterIterations;
         Outstates[1] = (double)primme.stats.numRestarts;
         Outstates[2] = (double)primme.stats.numMatvecs;
         Outstates[3] = (double)primme.stats.numPreconds; 
      }

      mxFree(evals);
      mxFree(evecs);
      mxFree(rnorms);
      primme_Free(&primme);
      if (outputfilename != NULL)
         mxFree(outputfilename);

   }
}

/******************************************************************************
 *
 * Matvec function. Callback to the appropriate function in MATLAB.
 *
 ******************************************************************************/

/* The MATLAB matvec handle is getMatvecHandle, which will receive the block
 * vector x, and perform y=A*x if A is a matrix, or y=afun(x) if the user 
 * provided the function afun. The resulting block vector y is returned to 
 * this function.
 */

void MatrixMatvec_d(void *x, void *y, int *blockSize, primme_params *primme)
{  
   double wt1 = primme_get_wtime(); 
   double * sendXr;
   double *xvec = (double *)x;
   mwSize n = primme->n;
   mwSize k;
   mwSize l;

   if (xvec == NULL) {
      mexErrMsgTxt("vector pointer x cannot be NULL pointer");
   }

   mxArray *rhs[1], *lhs[1];
   rhs[0] = mxCreateDoubleMatrix(n,*blockSize,mxREAL);
   sendXr = mxGetPr(rhs[0]);
   double * yvecr;
   double * ycopyvec = (double *)y;

   for (l = 0; l < n*(*blockSize); l++) {
      sendXr[l] = xvec[l];
   }

   mexCallMATLAB( 1, lhs, 1, rhs, "getMatvecHandle");        
   yvecr = mxGetPr(lhs[0]);

   for (l = 0; l < n*(*blockSize); l++) {
      ycopyvec[l] = yvecr[l];
   }

   mxDestroyArray(rhs[0]); 
   mxDestroyArray(lhs[0]);     

   double wt2 = primme_get_wtime(); 
   Matvec_mex_timer = Matvec_mex_timer + wt2 -wt1;

}


void MatrixMatvec_z(void *x, void *y, int *blockSize, primme_params *primme)
{  
   double wt1 = primme_get_wtime(); 

   double * sendXr;
   double * sendXi;
   Complex_Z *xvec = (Complex_Z *)x;
   mwSize n = primme->n;
   mwSize k;
   mwSize l;


   mxArray *rhs[1], *lhs[1];
   rhs[0] = mxCreateDoubleMatrix(n,*blockSize,mxCOMPLEX);
   sendXr = mxGetPr(rhs[0]);
   sendXi = mxGetPi(rhs[0]);
   double * yvecr;
   double * yveci;

   Complex_Z * ycopyvec = (Complex_Z *)y;

   for (l = 0; l < n*(*blockSize); l++) {
      sendXr[l] = xvec[l].r;
      sendXi[l] = xvec[l].i;
   }

   mexCallMATLAB( 1, lhs, 1, rhs, "getMatvecHandle");        
   yvecr = mxGetPr(lhs[0]);
   yveci = mxGetPi(lhs[0]);

   for (l = 0; l < n*(*blockSize); l++) {
      ycopyvec[l].r = yvecr[l];
      ycopyvec[l].i = yveci[l];
   }

   mxDestroyArray(rhs[0]);
   mxDestroyArray(lhs[0]);

   double wt2 = primme_get_wtime(); 
   Matvec_mex_timer = Matvec_mex_timer + wt2 -wt1;
}


/******************************************************************************
 *
 * Preconditioner function. Calls the necessary function call in MATLAB.
 *
 ******************************************************************************/

/* The MATLAB preconditioning handle is getPrecondHandle, which will receive 
 * the block * vector x, and perform y=P\x if P is a matrix, or y=afun(x) 
 * if the user provided the preconditioning function afun. The resulting block 
 * vector y is returned to this function.
 */

void Preconditioner_d(void *x, void *y, int *blockSize, primme_params *primme)
{ 

   double * sendXr;
   double *xvec = (double *)x;
   mwSize n = primme->n;
   mwSize k;
   mwSize l;

   if (xvec == NULL) {
      mexErrMsgTxt("vector pointer x cannot be NULL pointer");
   }

   mxArray *rhs[1], *lhs[1];
   rhs[0] = mxCreateDoubleMatrix(n,*blockSize,mxREAL);
   sendXr = mxGetPr(rhs[0]);
   double * yvecr;
   double * ycopyvec = (double *)y;

   for (l = 0; l < n*(*blockSize); l++) {
      sendXr[l] = xvec[l];
   }

   mexCallMATLAB( 1, lhs, 1, rhs, "getPrecondHandle");
   yvecr = mxGetPr(lhs[0]);

   for (l = 0; l < n*(*blockSize); l++) {
      ycopyvec[l] = yvecr[l];
   }

   mxDestroyArray(rhs[0]);
   mxDestroyArray(lhs[0]);

}

void Preconditioner_z(void *x, void *y, int *blockSize, primme_params *primme)
{ 
   double * sendXr;
   double * sendXi;
   Complex_Z *xvec = (Complex_Z *)x;
   mwSize n = primme->n;
   mwSize k;
   mwSize l;

   mxArray *rhs[1], *lhs[1];
   rhs[0] = mxCreateDoubleMatrix(n,*blockSize,mxCOMPLEX);
   sendXr = mxGetPr(rhs[0]);
   sendXi = mxGetPi(rhs[0]);
   double * yvecr;
   double * yveci;

   Complex_Z * ycopyvec = (Complex_Z *)y;

   for (l = 0; l < n*(*blockSize); l++) {
      sendXr[l] = xvec[l].r;
      sendXi[l] = xvec[l].i;
   }

   mexCallMATLAB( 1, lhs, 1, rhs, "getPrecondHandle");
   yvecr = mxGetPr(lhs[0]);
   yveci = mxGetPi(lhs[0]);

   for (l = 0; l < n*(*blockSize); l++) {
      ycopyvec[l].r = yvecr[l];
      ycopyvec[l].i = yveci[l];
   }

   mxDestroyArray(rhs[0]);
   mxDestroyArray(lhs[0]);
}


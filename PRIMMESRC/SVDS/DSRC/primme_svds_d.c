 
#include <stdlib.h>   /* mallocs, free */
#include <unistd.h>   /* gethostname */
#include <stdio.h>  
#include <math.h>  
#include "primme_svds.h"
#include "const.h"
#include "wtime.h"
#include "primme_svds_private_d.h"
#include "numerical_d.h"
 
/*******************************************************************************
 * Subroutine dprimme_svds - This routine is a front end used to perform 
 *    error checking on the input parameters, perform validation, 
 *    and make the call to main_iter. 
 *
 *    Calling dprimme_svds with all evals, evecs, resNorms set to NULL
 *    returns the int and real memory required in the following primme fields:
 *            int primme->intWorkSize : bytes of int workspace needed
 *       long int primme->realWorkSize: bytes of real workspace needed
 * 
 * INPUT/OUTPUT ARRAYS AND PARAMETERS
 * ----------------------------------
 * evals  Contains the converged Ritz values upon return.   Should be of size 
 *        primme->numEvals.
 * 
 * evecs  The local portions of the converged Ritz vectors.  The dimension of
 *        the array is at least primme->nLocal*primme->numEvals
 *
 * resNorms  The residual norms of the converged Ritz vectors.  Should be of 
 *           size primme->numEvals
 *  
 * primme  Structure containing various solver parameters and statistics
 *         See readme.txt for INPUT/OUTPUT variables in primme
 *
 * Return Value
 * ------------
 *  0 - Success
 * -1 - Failure to allocate workspace
 * -2 - Malloc failure in allocating a permutation integer array
 * -3 - main_iter encountered a problem
 * -4 ...-32 - Invalid input (parameters or primme struct) returned 
 *             by check_input()
 *
 ******************************************************************************/
 
int dprimme_svds(double *svals, double *svecs, double *resNorms, 
            primme_svds_params *primme_svds) {
      
   int ret;
   int i,j;
   int one = 1; /*used for Matvec*/
   int At = 0;
   const char notransp[] = "notransp"; /*used for Matvec*/
   const char transp[] = "transp"; /*used for Matvec*/
   double ztmp;
   double *realWork;
   primme_preset_method method;

   /* ----------------------------------------------------------- */
   /* Primme_svds_initialize must be called by users unless users */  
   /* specify all parameters in primme_svds structure. Check if   */
   /* primme_svds inputs are good for bounds, correct values etc. */
   /* ----------------------------------------------------------- */
   ret = primme_svds_check_input(svals, svecs, resNorms, primme_svds); 
   if (ret != 0) {
      fprintf(stderr,"primme_svds: check input failed - error code:%d\n",ret);
      return(ret);
   }

   /* -----------------------------------------------------------*/
   /* Set parameters for primme based on svds method when using  */
   /* normal equation method or hybrid method in the first stage */  
   /* -----------------------------------------------------------*/
   if (primme_svds->svdsMethod == primme_svds_normalequations ||
      primme_svds->svdsMethod == primme_svds_hybrid) {
      
      At = (primme_svds->n > primme_svds->m) ? 1 : 0;
      /* ----------------------------------------- */
      /* Set some defaults for sequential programs */
      /* ----------------------------------------- */
      if (primme_svds->numProcs == 1) {
         primme_svds->mLocal = primme_svds->m;
         primme_svds->nLocal = primme_svds->n;
         primme_svds->procID = 0;
         primme_svds->numProcs = 1;
      }

      /* -----------------------------------------------*/
      /* Initialize defaults in primme structure        */
      /* -----------------------------------------------*/
      primme_initialize(&primme_svds->primme);
  
      /* -----------------------------------------------*/
      /* Set important parameters for primme structure  */
      /* -----------------------------------------------*/
      primme_svds->primme.numEvals = primme_svds->numSvals;
      primme_svds->primme.aNorm = primme_svds->aNorm;
      primme_svds->primme.eps = primme_svds->eps;
      primme_svds->primme.correctionParams.precondition = 
                                      primme_svds->precondition;
      primme_svds->primme.initSize = primme_svds->initSize;
      primme_svds->primme.maxBasisSize = primme_svds->maxBasisSize;
      primme_svds->primme.maxBlockSize = primme_svds->maxBlockSize;
      primme_svds->primme.maxMatvecs = primme_svds->maxMatvecs;
      primme_svds->primme.iseed[0] = primme_svds->iseed[0];
      primme_svds->primme.iseed[1] = primme_svds->iseed[1];
      primme_svds->primme.iseed[2] = primme_svds->iseed[2];
      primme_svds->primme.iseed[3] = primme_svds->iseed[3];
      primme_svds->primme.printLevel = primme_svds->printLevel;
      primme_svds->primme.outputFile = primme_svds->outputFile;
      primme_svds->primme.matrix     = primme_svds;

      /* ---------------------------------------------- */
      /* Set some parameters only for parallel programs */
      /* ---------------------------------------------- */
      if (primme_svds->numProcs > 1 && primme_svds->globalSumDouble != NULL) {
           primme_svds->primme.procID = primme_svds->procID;
           primme_svds->primme.numProcs = primme_svds->numProcs;
           primme_svds->primme.commInfo = primme_svds->commInfo;
           primme_svds->primme.globalSumDouble = primme_svds->globalSumDouble;
      }

      primme_svds->primme.n = At == 0 ? primme_svds->n : primme_svds->m;
      primme_svds->primme.nLocal = At == 0 ? primme_svds->nLocal : primme_svds->mLocal;
      if (primme_svds->target == primme_svds_largest){
          primme_svds->primme.target = primme_largest;
      }
      else {
          primme_svds->primme.target = primme_smallest;  
      }
      method = primme_svds->eigsMethod_stage1;
      if (primme_set_method(method, &primme_svds->primme) < 0 ) {
         fprintf(primme_svds->primme.outputFile, 
                    "No preset method. Using custom settings\n");
      }
      primme_svds->primme.locking = 0;
      primme_svds->primme.projectionParams.projection = primme_RR;
      primme_svds->primme.projectionParams.refinedScheme = primme_DisableRef;
      if (primme_svds->primme.correctionParams.precondition){
          primme_svds->primme.InitBasisMode = 0;
      }
      else{
          primme_svds->primme.InitBasisMode = 1; /*lingfei: use 1 for debugging, later use 2 for general case*/
      }
      if(primme_svds->primme.aNorm > 0){
          if(primme_svds->svdsMethod == primme_svds_normalequations)
          {
              primme_svds->primme.aNorm = primme_svds->primme.aNorm*
                    primme_svds->primme.aNorm;
          }
      }
      if (primme_svds->svdsMethod == primme_svds_hybrid){
          primme_svds->primme.DefineConvCriteria = 1;//lingfei: may remove it later
          primme_svds->primme.ForceVerificationOnExit = 1;
      }

      /* ---------------------------------------------------*/
      /* Set up matrix vector and preconditioner for primme */
      /* -------------------------------------------------- */
      /* Allocate realWork array for MV and Precond operations */
      primme_svds->realWork = (double *)primme_calloc(
        primme_svds->numSvals*max(primme_svds->mLocal,primme_svds->nLocal), 
        sizeof(double), "primme_svds->realWork"); 
      if (primme_svds->realWork == NULL){
          fprintf(stderr,"primme_svds: realWork memory allocation failed\n");
          return MALLOC_FAILURE;
      }
      realWork = (double *)primme_svds->realWork;

      primme_svds->primme.matrixMatvec = MatrixATA_Matvec;
      primme_svds->primme.applyPreconditioner = MatrixATA_Precond;

      /* --------------------------------------------------------*/
      /* Call primme eigensolver for the first stage             */ 
      /* --------------------------------------------------------*/
      if (primme_svds->procID == 0)
         primme_display_params(primme_svds->primme);
     
      /* Normal equations or hybrid method is used, since we will 
         arrange U and V in order so when m>=n, use second part of
         memory for V; when m<n, use first part of memory for U.*/ 
      if (primme_svds->m >= primme_svds->n){
         ret = dprimme(svals, &svecs[primme_svds->numSvals*
            primme_svds->mLocal], resNorms, &primme_svds->primme); 
       }
       else{ 
          ret = dprimme(svals, svecs, resNorms, 
            &primme_svds->primme); 
       }
      if(ret != 0) {
        fprintf(stderr,"primme_svds: call PRIMME(ATA) failed - error code:%d\n",ret);
        return(CALL_PRIMME_ATA_FAILURE);
      }

      /* --------------------------------------------------------*/
      /* Record performance measurements from the first stage    */ 
      /* --------------------------------------------------------*/
      primme_svds->stats.numOuterIterations = 
            primme_svds->primme.stats.numOuterIterations;
      primme_svds->stats.numRestarts = primme_svds->primme.stats.numRestarts;
      primme_svds->stats.numMatvecs = primme_svds->primme.stats.numMatvecs;
      primme_svds->stats.numPreconds = primme_svds->primme.stats.numPreconds;
      primme_svds->stats.elapsedTime = primme_svds->primme.stats.elapsedTime;

      /* --------------------------------------------------------*/
      /* Formulate left and right singular vectors and perform:  */ 
      /* 1) if use normal equations method, return U, S, V       */
      /* 2) if use hybrid method in stage 1, pass the proper     */
      /*    initial guesses and shifts for hybrid in stage 2     */
      /* --------------------------------------------------------*/
      primme_svds->initSize = primme_svds->primme.initSize;
      if (primme_svds->m >= primme_svds->n){
          for (i=0;i<primme_svds->initSize;i++){
              svals[i] = sqrt(fabs(svals[i]));
              (* primme_svds->matrixMatvec)(&svecs[
                 primme_svds->initSize*primme_svds->mLocal + i*primme_svds->nLocal],
                  &svecs[i*primme_svds->mLocal], &one, primme_svds,notransp);
                  ztmp = 1.0L/svals[i];
                  Num_scal_dprimme(primme_svds->mLocal, ztmp,
                    &svecs[i*primme_svds->mLocal], 1);
          }
      }
      else {
          for (i=0;i<primme_svds->initSize;i++){
              svals[i] = sqrt(fabs(svals[i]));
              (* primme_svds->matrixMatvec)(&svecs
                [i*primme_svds->mLocal], &svecs[primme_svds->initSize*
                primme_svds->mLocal + i*primme_svds->nLocal], &one, primme_svds,transp);
                  ztmp = 1.0L/svals[i];
                  Num_scal_dprimme(primme_svds->nLocal, ztmp, 
                            &svecs[primme_svds->initSize*primme_svds->mLocal 
                            + i*primme_svds->nLocal], 1);
           }
      }
      primme_svds->aNorm = sqrt(primme_svds->primme.aNorm);
      primme_Free(&primme_svds->primme);
      if (primme_svds->svdsMethod == primme_svds_normalequations){
          return 0;
      }
      primme_svds_Free(primme_svds);
   }
   
   /* ----------------------------------------------------------*/
   /* Set parameters for primme based on svds method when using */
   /* augmented method or hybrid method in the second stage     */  
   /* ----------------------------------------------------------*/
   if (primme_svds->svdsMethod == primme_svds_augmented || 
       primme_svds->svdsMethod == primme_svds_hybrid) {
      /* ----------------------------------------- */
      /* Set some defaults for sequential programs */
      /* ----------------------------------------- */
      if (primme_svds->numProcs == 1) {
         primme_svds->mLocal = primme_svds->m; 
         primme_svds->nLocal = primme_svds->n; 
         primme_svds->procID = 0;
         primme_svds->numProcs = 1;
      }

      /* -----------------------------------------------*/
      /* Initialize defaults in primme structure        */
      /* -----------------------------------------------*/
      primme_initialize(&primme_svds->primme);
  
      /* -----------------------------------------------*/
      /* Set important parameters for primme structure  */
      /* -----------------------------------------------*/
      primme_svds->primme.numEvals = primme_svds->numSvals;
      primme_svds->primme.aNorm = primme_svds->aNorm;
      primme_svds->primme.eps = primme_svds->eps;
      primme_svds->primme.correctionParams.precondition = 
                                      primme_svds->precondition;
      primme_svds->primme.initSize = primme_svds->initSize;
      primme_svds->primme.maxBasisSize = primme_svds->maxBasisSize;
      if (primme_svds->target == primme_svds_largest)
          primme_svds->primme.maxBlockSize = primme_svds->maxBlockSize;
      else { /*primme_svds->target == primme_svds_smallest*/
          if (primme_svds->maxBlockSize > 1){
              primme_svds->maxBlockSize = 1; /*Only support block size 1*/
              primme_svds->primme.maxBlockSize = 1;
          }
      }
      primme_svds->primme.maxMatvecs = primme_svds->maxMatvecs;
      primme_svds->primme.iseed[0] = primme_svds->iseed[0];
      primme_svds->primme.iseed[1] = primme_svds->iseed[1];
      primme_svds->primme.iseed[2] = primme_svds->iseed[2];
      primme_svds->primme.iseed[3] = primme_svds->iseed[3];
      primme_svds->primme.printLevel = primme_svds->printLevel;
      primme_svds->primme.outputFile = primme_svds->outputFile;
      primme_svds->primme.matrix     = primme_svds;

      /* ---------------------------------------------- */
      /* Set some parameters only for parallel programs */
      /* ---------------------------------------------- */
      if (primme_svds->numProcs > 1 && primme_svds->globalSumDouble != NULL) {
           primme_svds->primme.procID = primme_svds->procID;
           primme_svds->primme.numProcs = primme_svds->numProcs;
           primme_svds->primme.commInfo = primme_svds->commInfo;
           primme_svds->primme.globalSumDouble = primme_svds->globalSumDouble;
      }
      
      primme_svds->primme.n = primme_svds->m+primme_svds->n;
      primme_svds->primme.nLocal = primme_svds->mLocal+primme_svds->nLocal;
      if (primme_svds->target == primme_svds_largest) {
          primme_svds->primme.target = primme_largest;
          method = primme_svds->eigsMethod_stage2;
          if (primme_set_method(method, &primme_svds->primme) < 0 ) {
              fprintf(primme_svds->primme.outputFile, 
                    "No preset method. Using custom settings\n");
          }
          primme_svds->primme.projectionParams.projection = primme_RR;
          primme_svds->primme.projectionParams.refinedScheme = primme_DisableRef;
          primme_svds->primme.locking = 0;
          if (primme_svds->svdsMethod == primme_svds_hybrid){
             primme_svds->primme.InitBasisMode = 1;
//             primme_svds->primme.initSize = primme_svds->numSvals;
          }
      }
      else {
          primme_svds->primme.target = primme_closest_geq;  
          method = primme_svds->eigsMethod_stage2;
          if (primme_set_method(method, &primme_svds->primme) < 0 ) {
              fprintf(primme_svds->primme.outputFile, 
                    "No preset method. Using custom settings\n");
          }
          primme_svds->primme.locking = 1;
          primme_svds->primme.projectionParams.projection = primme_RR_Refined;
          primme_svds->primme.projectionParams.refinedScheme = primme_OneAccuShift_QR;
          if (primme_svds->svdsMethod == primme_svds_augmented){
              primme_svds->primme.numTargetShifts = primme_svds->numTargetShifts;
              if (primme_svds->primme.numTargetShifts > 0) {
                  primme_svds->primme.targetShifts = (double *)primme_calloc(
                        primme_svds->primme.numTargetShifts, sizeof(double), 
                        "targetShifts");
                  for (i=0; i<primme_svds->primme.numTargetShifts;i++) {
                      primme_svds->primme.targetShifts[i] = 
                                                primme_svds->targetShifts[i];
                  }
              }
              primme_svds->primme.qr_need = 1;
              primme_svds->primme.InitBasisMode = 1;
              primme_svds->primme.ReIntroInitGuessToBasis = 0;
          }
          else if (primme_svds->svdsMethod == primme_svds_hybrid){
              primme_svds->primme.numTargetShifts = primme_svds->numSvals;
              primme_svds->primme.targetShifts = (double *)primme_calloc(
                    primme_svds->primme.numTargetShifts, sizeof(double), 
                    "targetShifts");
              for(i=0; i< primme_svds->primme.numTargetShifts; i++){
                  primme_svds->primme.targetShifts[i] = svals[i];
              }
              primme_svds->primme.qr_need = 1;
              primme_svds->primme.InitBasisMode = 1;
              primme_svds->primme.ReIntroInitGuessToBasis = 1;
//              primme_svds->primme.initSize = primme_svds->numSvals;
          }
      }
      /* ---------------------------------------------------*/
      /* Set up matrix vector and preconditioner for primme */
      /* -------------------------------------------------- */
      /* Allocate realWork array for MV and Precond operations */
      primme_svds->realWork = (double *)primme_calloc(
        primme_svds->numSvals*(primme_svds->mLocal+primme_svds->nLocal), 
        sizeof(double), "primme_svds->realWork"); 
      if (primme_svds->realWork == NULL){
          fprintf(stderr,"primme_svds: realWork memory allocation failed\n");
          return MALLOC_FAILURE;
      }
      realWork = (double *)primme_svds->realWork;

      primme_svds->primme.matrixMatvec = MatrixB_Matvec;
      primme_svds->primme.applyPreconditioner = MatrixB_Precond;

      /* ---------------------------------------------------*/
      /* Set up initial vectors for PRIMME                  */
      /* -------------------------------------------------- */
      int num = min(primme_svds->initSize, primme_svds->numSvals);
      Num_dcopy_dprimme(num*primme_svds->mLocal, svecs, 1, realWork, 1);
      Num_dcopy_dprimme(num*primme_svds->nLocal, &svecs[primme_svds->initSize*primme_svds->mLocal],
                                  1, &realWork[num*primme_svds->mLocal], 1);
      for (i=0; i<num; i++) {
          Num_dcopy_dprimme(primme_svds->nLocal,
                &realWork[num*primme_svds->mLocal+i*primme_svds->nLocal], 1,
                &svecs[i*(primme_svds->nLocal+primme_svds->mLocal)], 1);
      }
      for (i=0; i<num; i++) {
          Num_dcopy_dprimme(primme_svds->mLocal, 
                &realWork[i*primme_svds->mLocal], 1,
                &svecs[i*(primme_svds->nLocal+primme_svds->mLocal) 
                + primme_svds->nLocal], 1);
      }
      
      /* --------------------------------------------------------*/
      /* Call primme eigensolver for the second stage            */ 
      /* --------------------------------------------------------*/
      if (primme_svds->procID == 0)
         primme_display_params(primme_svds->primme);
     
      ret = dprimme(svals, svecs, resNorms, &primme_svds->primme); 
      
      if(ret != 0) {
        fprintf(stderr,"primme_svds: call PRIMME(B) failed - error code:%d\n",ret);
        return(CALL_PRIMME_B_FAILURE);
      }
      
      /* --------------------------------------------------------*/
      /* Record performance measurements from the second stage   */ 
      /* --------------------------------------------------------*/
      primme_svds->stats.numOuterIterations += 
            primme_svds->primme.stats.numOuterIterations;
      primme_svds->stats.numRestarts += primme_svds->primme.stats.numRestarts;
      primme_svds->stats.numMatvecs += primme_svds->primme.stats.numMatvecs;
      primme_svds->stats.numPreconds += primme_svds->primme.stats.numPreconds;
      primme_svds->stats.elapsedTime += primme_svds->primme.stats.elapsedTime;

      /* -------------------------------------------------------------*/
      /* Scale svecs by sqrt(2) to obtain the left and right singular */
      /* vectors and formulate U and V in order                       */ 
      /* -------------------------------------------------------------*/
      ztmp = sqrt(2.0L);
      Num_scal_dprimme((primme_svds->mLocal + primme_svds->nLocal)*
            primme_svds->numSvals, ztmp, svecs, 1);

      /* Suppose returned eigenvectors are [v1 u1 v2 u2 ... vk uk],
         so we copy v1, v2, ..., vk to primme_svds->realWork, and 
         move u1, u2, ..., uk to the first part of svecs memory. 
         Then we move primme_svds->realWork to the second part of 
         svecs memory for V.  */
      for (i=0;i<primme_svds->numSvals;i++){
          Num_dcopy_dprimme(primme_svds->nLocal, &svecs[i*
                (primme_svds->nLocal + primme_svds->mLocal)], 1, &realWork[i*
                primme_svds->nLocal], 1); 
          Num_dcopy_dprimme(primme_svds->mLocal, &svecs[i*
                (primme_svds->mLocal+primme_svds->nLocal) + primme_svds->nLocal], 1, 
                &svecs[i*(primme_svds->nLocal+ primme_svds->mLocal)], 1); 
      }
      Num_dcopy_dprimme(primme_svds->nLocal*primme_svds->numSvals, 
            realWork, 1, &svecs[primme_svds->numSvals*primme_svds->mLocal], 1); 
      primme_svds->aNorm = primme_svds->primme.aNorm;
      if (primme_svds->procID == 0)
         primme_display_params(primme_svds->primme);//lingfei: remove it later
      primme_Free(&primme_svds->primme);
   }

   return(0);
}

/******************************************************************************
 *
 * static int primme_svds_check_input(double *svals, double *svecs, double *resNorms, 
 *                        primme_svds_params *primme_svds) 
 *
 * INPUT
 * -----
 *  svals, svecs, resNorms   Output arrays for primme
 *  primme_svds              the main structure of parameters 
 *
 * return value -   0    If input parameters in primme are appropriate
 *              -4..-19  Inappropriate input parameters were found
 *
 ******************************************************************************/
static int primme_svds_check_input(double *svals, double *svecs, double *resNorms, 
                       primme_svds_params *primme_svds) {
   int ret;
   ret = 0;

   if (primme_svds == NULL)
      ret = -4;
   else if (primme_svds->n <= 0 || primme_svds->m <= 0) 
      ret = -5;
   else if (primme_svds->numProcs < 1)
      ret = -6;
   else if (primme_svds->matrixMatvec == NULL) 
      ret = -7;
   else if (primme_svds->applyPreconditioner == NULL && 
       primme_svds->precondition != 0) 
      ret = -8;
   else if (primme_svds->matrix == NULL) 
      ret = -9;
   else if (primme_svds->numProcs >1 && primme_svds->globalSumDouble == NULL)
      ret = -10;
   else if (primme_svds->numSvals > primme_svds->n)
      ret = -11;
   else if (primme_svds->numSvals < 1)
      ret = -12;
   else if (primme_svds->maxBlockSize < 1)
      ret = -13;
   else if ( primme_svds->target != primme_svds_smallest  &&
             primme_svds->target != primme_svds_largest)
      ret = -14;
   else if ( primme_svds->svdsMethod != primme_svds_hybrid &&
             primme_svds->svdsMethod != primme_svds_normalequations &&
             primme_svds->svdsMethod != primme_svds_augmented)
      ret = -15;
   else if (primme_svds->printLevel < 0 || primme_svds->printLevel > 5)
      ret = -16; 
   else if (svals == NULL)
      ret = -17;
   else if (svecs == NULL)
      ret = -18;
   else if (resNorms == NULL)
      ret = -19;

   return ret;
  /***************************************************************************/
} /* end of check_input
   ***************************************************************************/

/**********************************************************************************
 * void MatrixATA_Matvec(void *x, void *y, int *blockSize, primme_params *primme) *
 **********************************************************************************/
void MatrixATA_Matvec(void *x, void *y, int *blockSize, primme_params *primme){
    /*When primme_sdvs->m >= primme_svds->n, perform A'(Ax) */
    const char notransp[] = "notransp"; /*used for Matvec*/
    const char transp[] = "transp"; /*used for Matvec*/
    primme_svds_params *primme_svds;
    primme_svds = (primme_svds_params *) primme->matrix;

    int one = 1; 
    int i; 
    double *xvec;
    double *yvec;
    void *xcopy;
    void *ycopy;
    
    xvec = (double *)x;
    yvec = (double *)y;
       
    /* To better support blocked Matvec operation for square and rectangular 
       matrix, it is better to call user Matvec function blockSize rather
       than calling user Matvec function once with a for loop. It greatly
       reduce the possibilities for all kinds of memory issues. */
    if(primme_svds->m >= primme_svds->n){
        for (i=0;i<*blockSize;i++) {
            xcopy = &xvec[i*primme->nLocal];
            ycopy = &yvec[i*primme->nLocal];
            (* primme_svds->matrixMatvec)(xcopy, primme_svds->realWork, &one, 
                                                    primme_svds, notransp);
            (* primme_svds->matrixMatvec)(primme_svds->realWork, ycopy, &one, 
                                                    primme_svds, transp);
        }
    }
    else{ /* primme_svds->m < primme_svds->n */
        for (i=0;i<*blockSize;i++) {
            xcopy = &xvec[i*primme->nLocal];
            ycopy = &yvec[i*primme->nLocal];
            (* primme_svds->matrixMatvec)(xcopy, primme_svds->realWork, &one, 
                                                    primme_svds, transp);
            (* primme_svds->matrixMatvec)(primme_svds->realWork, ycopy, &one, 
                                                    primme_svds, notransp);
        }
    }
}

/**********************************************************************************
 * void MatrixB_Matvec(void *x, void *y, int *blockSize, primme_params *primme) *
 **********************************************************************************/
void MatrixB_Matvec(void *x, void *y, int *blockSize, primme_params *primme){
    /*When performing [0 A'; A 0]x, perform A'x(1+N:M+N,:) and Ax(1:N,:) */
    /*lingfei: when blockSize is not 1, we have trouble. Resolve it later*/
    const char notransp[] = "notransp"; /*used for Matvec*/
    const char transp[] = "transp"; /*used for Matvec*/
    primme_svds_params *primme_svds;
    primme_svds = (primme_svds_params *) primme->matrix;
    double *xvec;
    double *yvec;
    void *xcopy;
    void *ycopy;

    xvec = (double *)x;
    yvec = (double *)y;
    xcopy = &xvec[primme_svds->nLocal];
    ycopy = &yvec[primme_svds->nLocal];
    
    if((primme_svds->target == primme_svds_smallest && 
        primme_svds->maxBlockSize == 1) || 
        primme_svds->target == primme_svds_largest){
        (* primme_svds->matrixMatvec)(xcopy, y, blockSize, primme_svds, transp);
        (* primme_svds->matrixMatvec)(x, ycopy, blockSize, primme_svds, notransp);
    }
    else
      fprintf(stderr,"primme_svds: maxBlockSize must be 1 when using primme_svds_smallest\n");

}

/**********************************************************************************
 * void MatrixATA_Precond(void *x, void *y, int *blockSize, primme_params *primme) *
 **********************************************************************************/
void MatrixATA_Precond(void *x, void *y, int *blockSize, primme_params *primme){
    /*When primme_sdvs->m >= primme_svds->n, perform A'(Ax) */
   const char notransp[] = "notransp"; /*used for Matvec*/
   const char transp[] = "transp"; /*used for Matvec*/
    primme_svds_params *primme_svds;
    primme_svds = (primme_svds_params *) primme->matrix;
        
    int one = 1; 
    int i; 
    double *xvec;
    double *yvec;
    void *xcopy;
    void *ycopy;
    
    xvec = (double *)x;
    yvec = (double *)y;
       
    /* To better support blocked preconditioning operation for square and rectangular 
       matrix, it is better to call user preconditioning function blockSize rather
       than calling user preconditioning function once with a for loop. It greatly
       reduce the possibilities for all kinds of memory issues. */
    if(primme_svds->m >= primme_svds->n){
        for (i=0;i<*blockSize;i++) {
            xcopy = &xvec[i*primme->nLocal];
            ycopy = &yvec[i*primme->nLocal];
            (* primme_svds->applyPreconditioner)(xcopy, primme_svds->realWork, &one, 
                                                    primme_svds, transp);
            (* primme_svds->applyPreconditioner)(primme_svds->realWork, ycopy, &one, 
                                                    primme_svds, notransp);
        }
    }
    else{ /* primme_svds->m < primme_svds->n */
        for (i=0;i<*blockSize;i++) {
            xcopy = &xvec[i*primme->nLocal];
            ycopy = &yvec[i*primme->nLocal];
            (* primme_svds->applyPreconditioner)(xcopy, primme_svds->realWork, &one, 
                                                    primme_svds, notransp);
            (* primme_svds->applyPreconditioner)(primme_svds->realWork, ycopy, &one, 
                                                    primme_svds, transp);
        }
    }
}

/**********************************************************************************
 * void MatrixB_Precond(void *x, void *y, int *blockSize, primme_params *primme) *
 **********************************************************************************/
void MatrixB_Precond(void *x, void *y, int *blockSize, primme_params *primme){
    /*When performing [0 A'; A 0]x, perform A'x(1+N:M+N,:) and Ax(1:N,:) */
    /*lingfei: when blockSize is not 1, we have trouble. Resolve it later*/
    const char notransp[] = "notransp"; /*used for Matvec*/
    const char transp[] = "transp"; /*used for Matvec*/
    primme_svds_params *primme_svds;
    primme_svds = (primme_svds_params *) primme->matrix;
    double *xvec;
    double *yvec;
    void *xcopy;
    void *ycopy;

    xvec = (double *)x;
    yvec = (double *)y;
    xcopy = &xvec[primme_svds->nLocal];
    ycopy = &yvec[primme_svds->nLocal];

    if((primme_svds->target == primme_svds_smallest && 
        primme_svds->maxBlockSize == 1) || 
        primme_svds->target == primme_svds_largest){
        (* primme_svds->applyPreconditioner)(xcopy, y, blockSize, 
                                                    primme_svds, notransp);
        (* primme_svds->applyPreconditioner)(x, ycopy, blockSize, 
                                                    primme_svds, transp);
    }
    else
      fprintf(stderr,"primme_svds: maxBlockSize must be 1 when using primme_svds_smallest\n");
    
}

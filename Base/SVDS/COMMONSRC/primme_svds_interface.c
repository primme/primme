
#include <stdlib.h>   /* mallocs, free */
#include <unistd.h>   /* gethostname */
#include <stdio.h>    
#include "primme_svds.h"

/***************************************************************************

   Initialize the primme_svds data structure
  
***************************************************************************/
void primme_svds_initialize(primme_svds_params *primme_svds) {

   /* Essential parameters */
   primme_svds->n                       = 0;
   primme_svds->m                       = 0;
   primme_svds->numSvals                = 6;
   primme_svds->target                  = primme_svds_largest;
   primme_svds->eigsMethod_stage1       = DEFAULT_MIN_MATVECS;
   primme_svds->eigsMethod_stage2       = JDQMR;
   primme_svds->svdsMethod              = primme_svds_hybrid;

   /* Shifts for primme_svds_augmented method */
   primme_svds->numTargetShifts         = 0;
   primme_svds->targetShifts            = NULL;

   /* Parallel computing parameters */
   primme_svds->numProcs                = 1;
   primme_svds->procID                  = 0;
   primme_svds->nLocal                  = 0;
   primme_svds->commInfo                = NULL;
   primme_svds->globalSumDouble         = NULL;

   /* Use these pointers for d/zprimme_svds function */
   primme_svds->primme_svds_info        = NULL;
   primme_svds->realWork                = NULL;

   /* Use these pointers to provide matrix/preconditioner */
   primme_svds->matrix                  = NULL;
   primme_svds->preconditioner          = NULL;

   /* Matvec and preconditioner */
   primme_svds->matrixMatvec            = NULL;
   primme_svds->applyPreconditioner     = NULL;

   /* Other important parameters users may set */
   primme_svds->aNorm                   = 0.0L;
   primme_svds->eps                     = 1e-12;
   primme_svds->precondition            = 0;
   primme_svds->initSize                = 0;
   primme_svds->maxBasisSize            = 0;
   primme_svds->maxBlockSize            = 1;
   primme_svds->maxMatvecs              = INT_MAX;
   primme_svds->printLevel              = 1;
   primme_svds->outputFile              = stdout;

   /* Reporting performance */
   primme_svds->stats.numOuterIterations= 0;
   primme_svds->stats.numRestarts       = 0;
   primme_svds->stats.numMatvecs        = 0;
   primme_svds->stats.numPreconds       = 0;
   primme_svds->stats.elapsedTime       = 0.0L;

   /* Internally used variables */
   primme_svds->iseed[0] = -1;   /* To set iseed, we first need procID           */                                                
   primme_svds->iseed[1] = -1;   /* Thus we set all iseeds to -1                 */
   primme_svds->iseed[2] = -1;   /* Unless users provide their own iseeds,       */
   primme_svds->iseed[3] = -1;   /* PRIMME will set thse later uniquely per proc */

}

/******************************************************************************
 *
 * void primme_svds_display_params(primme_params *primme);
 *
 *    Displays the current configuration of primme data structure
 *
 *****************************************************************************/
void primme_svds_display_params(primme_svds_params primme_svds) {

int i;
FILE *outputFile = primme_svds.outputFile;

fprintf(outputFile, "// ---------------------------------------------------\n");
fprintf(outputFile, "//            primme_svds configuration               \n");
fprintf(outputFile, "// ---------------------------------------------------\n");

fprintf(outputFile, "primme_svds.m = %d \n",primme_svds.m);
fprintf(outputFile, "primme_svds.n = %d \n",primme_svds.n);
fprintf(outputFile, "primme_svds.nLocal = %d \n",primme_svds.nLocal);
fprintf(outputFile, "primme_svds.numProcs = %d \n",primme_svds.numProcs);
fprintf(outputFile, "primme_svds.procID = %d \n",primme_svds.procID);

fprintf(outputFile, "\n// Output and reporting\n");
fprintf(outputFile, "primme_svds.printLevel = %d \n",primme_svds.printLevel);

fprintf(outputFile, "\n// Solver parameters\n");
fprintf(outputFile, "primme_svds.numSvals = %d \n",primme_svds.numSvals);
switch (primme_svds.target){
   case primme_svds_smallest:
      fprintf(outputFile, "primme_svds.target = primme_svds_smallest\n");
      break;
   case primme_svds_largest:
      fprintf(outputFile, "primme_svds.target = primme_svds_largest\n");
      break;
}
switch (primme_svds.svdsMethod){
   case primme_svds_hybrid:
      fprintf(outputFile, "primme_svds.svdsMethod = primme_svds_hybrid\n");
      break;
   case primme_svds_normalequations:
      fprintf(outputFile, "primme_svds.svdsMethod = primme_svds_normalequations\n");
      break;
   case primme_svds_augmented:
      fprintf(outputFile, "primme_svds.svdsMethod = primme_svds_augmented\n");
      break;
}
switch (primme_svds.eigsMethod_stage1){
   case DYNAMIC: 
      fprintf(outputFile, "primme_svds.eigsMethod_stage2 = DYNAMIC\n");
      break;
   case DEFAULT_MIN_TIME:
      fprintf(outputFile, "primme_svds.eigsMethod_stage1 = DEFAULT_MIN_TIME\n");
      break;
   case DEFAULT_MIN_MATVECS:
      fprintf(outputFile, "primme_svds.eigsMethod_stage1 = DEFAULT_MIN_MATVECS\n");
      break;
   case Arnoldi:
      fprintf(outputFile, "primme_svds.eigsMethod_stage1 = Arnoldi\n");
      break;
   case GD:
      fprintf(outputFile, "primme_svds.eigsMethod_stage1 = GD\n");
      break;
   case GD_plusK:
      fprintf(outputFile, "primme_svds.eigsMethod_stage1 = GD_plusK\n");
      break;
   case GD_Olsen_plusK:
      fprintf(outputFile, "primme_svds.eigsMethod_stage1 = GD_Olsen_plusK\n");
      break;
   case JD_Olsen_plusK:
      fprintf(outputFile, "primme_svds.eigsMethod_stage1 = JD_Olsen_plusK\n");
      break;
   case RQI:
      fprintf(outputFile, "primme_svds.eigsMethod_stage1 = RQI\n");
      break;
   case JDQR:
      fprintf(outputFile, "primme_svds.eigsMethod_stage1 = JDQR\n");
      break;
   case JDQMR:
      fprintf(outputFile, "primme_svds.eigsMethod_stage1 = JDQMR\n");
      break;
   case JDQMR_ETol:
      fprintf(outputFile, "primme_svds.eigsMethod_stage1 = JDQMR_ETol\n");
      break;
   case SUBSPACE_ITERATION:
      fprintf(outputFile, "primme_svds.eigsMethod_stage1 = SUBSPACE_ITERATION\n");
      break;
   case LOBPCG_OrthoBasis:
      fprintf(outputFile, "primme_svds.eigsMethod_stage1 = LOBPCG_OrthoBasis\n");
      break;
   case LOBPCG_OrthoBasis_Window:
      fprintf(outputFile, "primme_svds.eigsMethod_stage1 = LOBPCG_OrthoBasis_Window\n");
      break;
}
switch (primme_svds.eigsMethod_stage2){
   case DYNAMIC: 
      fprintf(outputFile, "primme_svds.eigsMethod_stage2 = DYNAMIC\n");
      break;
   case DEFAULT_MIN_TIME:
      fprintf(outputFile, "primme_svds.eigsMethod_stage2 = DEFAULT_MIN_TIME\n");
      break;
   case DEFAULT_MIN_MATVECS:
      fprintf(outputFile, "primme_svds.eigsMethod_stage2 = DEFAULT_MIN_MATVECS\n");
      break;
   case Arnoldi:
      fprintf(outputFile, "primme_svds.eigsMethod_stage2 = Arnoldi\n");
      break;
   case GD:
      fprintf(outputFile, "primme_svds.eigsMethod_stage2 = GD\n");
      break;
   case GD_plusK:
      fprintf(outputFile, "primme_svds.eigsMethod_stage2 = GD_plusK\n");
      break;
   case GD_Olsen_plusK:
      fprintf(outputFile, "primme_svds.eigsMethod_stage2 = GD_Olsen_plusK\n");
      break;
   case JD_Olsen_plusK:
      fprintf(outputFile, "primme_svds.eigsMethod_stage2 = JD_Olsen_plusK\n");
      break;
   case RQI:
      fprintf(outputFile, "primme_svds.eigsMethod_stage2 = RQI\n");
      break;
   case JDQR:
      fprintf(outputFile, "primme_svds.eigsMethod_stage2 = JDQR\n");
      break;
   case JDQMR:
      fprintf(outputFile, "primme_svds.eigsMethod_stage2 = JDQMR\n");
      break;
   case JDQMR_ETol:
      fprintf(outputFile, "primme_svds.eigsMethod_stage2 = JDQMR_ETol\n");
      break;
   case SUBSPACE_ITERATION:
      fprintf(outputFile, "primme_svds.eigsMethod_stage2 = SUBSPACE_ITERATION\n");
      break;
   case LOBPCG_OrthoBasis:
      fprintf(outputFile, "primme_svds.eigsMethod_stage2 = LOBPCG_OrthoBasis\n");
      break;
   case LOBPCG_OrthoBasis_Window:
      fprintf(outputFile, "primme_svds.eigsMethod_stage2 = LOBPCG_OrthoBasis_Window\n");
      break;
}
fprintf(outputFile, "primme_svds.aNorm = %e \n",primme_svds.aNorm);
fprintf(outputFile, "primme_svds.eps = %e \n",primme_svds.eps);
fprintf(outputFile, "primme_svds.maxBasisSize = %d \n",primme_svds.maxBasisSize);
fprintf(outputFile, "primme_svds.maxBlockSize = %d\n",primme_svds.maxBlockSize);
fprintf(outputFile, "primme_svds.maxMatvecs = %d\n",primme_svds.maxMatvecs);
fprintf(outputFile, "primme_svds.initSize = %d\n",primme_svds.initSize);
fprintf(outputFile, "primme_svds.precondition = %d\n",primme_svds.precondition);
fprintf(outputFile, "primme_svds.iseed =");
for (i=0; i<4;i++) {
    fprintf(outputFile, " %d",primme_svds.iseed[i]);
}
fprintf(outputFile, "\n");
fprintf(outputFile, "primme_svds.numTargetShifts = %d\n",primme_svds.numTargetShifts);    
if (primme_svds.numTargetShifts > 0) {
    fprintf(outputFile, "primme_svds.targetShifts =");
    for (i=0; i<primme_svds.numTargetShifts;i++) {
        fprintf(outputFile, " %e",primme_svds.targetShifts[i]);
    }   
    fprintf(outputFile, "\n");
}

fprintf(outputFile, "\n");
fprintf(outputFile, "// ---------------------------------------------------\n");
fflush(outputFile);

  /**************************************************************************/
} /* end of display params */


void primme_svds_Free(primme_svds_params *params) {
    
    free(params->realWork);  
}


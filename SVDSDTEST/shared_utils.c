#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include "shared_utils.h"

/******************************************************************************
 *
 * Reads the solver parameters for configuring the primme_svds structure of 
 * dprimme_svds(). 
 *
******************************************************************************/
int read_solver_params(char *configFileName, char *outputFileName,
		primme_svds_params *primme_svds) {

   int line, ret, i;
   char ident[128];
   char op[128];
   char stringValue[128];
   FILE *configFile;

   if ((configFile = fopen(configFileName, "r")) == NULL) {
      fprintf(stderr,"Error(read_solver_params): Could not open config file\n");
      fprintf(stderr, "config file: %s\n", configFileName);
      return(-1);
   }
    

   line = 1;
   while (EOF != fscanf(configFile, "%s", ident)) {
      if (strcmp(ident, "//") == 0) {
         fgets(ident, 2048, configFile);
         line++;
         continue;
      }
      else {
         ret = fscanf(configFile, "%s", op);
      }

      if (ret == EOF) {
         fprintf(stderr, "ERROR(read_solver_params): Unexpected end of file\n");
         return(-1);
      }

      if (strcmp(op, "=") == 0) {
         if (strcmp(ident, "primme_svds.eigsMethod_stage1") == 0) { 
            ret = fscanf(configFile, "%s\n", stringValue);
            if (ret == 1) {
               if (strcmp(stringValue,      "DYNAMIC") == 0) {
		       primme_svds->eigsMethod_stage1 = DYNAMIC;
               }
	           else if (strcmp(stringValue, "DEFAULT_MIN_TIME") == 0) {
		       primme_svds->eigsMethod_stage1 = DEFAULT_MIN_TIME;
               }
               else if (strcmp(stringValue, "DEFAULT_MIN_MATVECS") == 0) {
		       primme_svds->eigsMethod_stage1 = DEFAULT_MIN_MATVECS;
               }
               else if (strcmp(stringValue, "Arnoldi") == 0) {
		       primme_svds->eigsMethod_stage1 = Arnoldi;
               }
               else if (strcmp(stringValue, "GD") == 0) {
		       primme_svds->eigsMethod_stage1 = GD;
               }
               else if (strcmp(stringValue, "GD_plusK") == 0) {
		       primme_svds->eigsMethod_stage1 = GD_plusK;
               }
               else if (strcmp(stringValue, "GD_Olsen_plusK") == 0) {
		       primme_svds->eigsMethod_stage1 = GD_Olsen_plusK;
               }
               else if (strcmp(stringValue, "JD_Olsen_plusK") == 0) {
		       primme_svds->eigsMethod_stage1 = JD_Olsen_plusK;
               }
               else if (strcmp(stringValue, "RQI") == 0) {
		       primme_svds->eigsMethod_stage1 = RQI;
               }
               else if (strcmp(stringValue, "JDQR") == 0) {
		       primme_svds->eigsMethod_stage1 = JDQR;
               }
               else if (strcmp(stringValue, "JDQMR") == 0) {
		       primme_svds->eigsMethod_stage1 = JDQMR;
               }
               else if (strcmp(stringValue, "JDQMR_ETol") == 0) {
		       primme_svds->eigsMethod_stage1 = JDQMR_ETol;
               }
               else if (strcmp(stringValue, "SUBSPACE_ITERATION") == 0) {
		       primme_svds->eigsMethod_stage1 = SUBSPACE_ITERATION;
               }
               else if (strcmp(stringValue, "LOBPCG_OrthoBasis") == 0) {
		       primme_svds->eigsMethod_stage1 = LOBPCG_OrthoBasis;
               }
               else if (strcmp(stringValue, "LOBPCG_OrthoBasis_Window") == 0) {
		       primme_svds->eigsMethod_stage1 = LOBPCG_OrthoBasis_Window;
               }
               else {
                  printf("Invalid eigenMethod value\n");
                  ret = 0;
               }
            }
         }
         else if (strcmp(ident, "primme_svds.eigsMethod_stage2") == 0){
            ret = fscanf(configFile, "%s\n", stringValue);
            if (ret == 1) {
               if (strcmp(stringValue,      "DYNAMIC") == 0) {
		       primme_svds->eigsMethod_stage2 = DYNAMIC;
               }
	           else if (strcmp(stringValue, "DEFAULT_MIN_TIME") == 0) {
		       primme_svds->eigsMethod_stage2 = DEFAULT_MIN_TIME;
               }
               else if (strcmp(stringValue, "DEFAULT_MIN_MATVECS") == 0) {
		       primme_svds->eigsMethod_stage2 = DEFAULT_MIN_MATVECS;
               }
               else if (strcmp(stringValue, "Arnoldi") == 0) {
		       primme_svds->eigsMethod_stage2 = Arnoldi;
               }
               else if (strcmp(stringValue, "GD") == 0) {
		       primme_svds->eigsMethod_stage2 = GD;
               }
               else if (strcmp(stringValue, "GD_plusK") == 0) {
		       primme_svds->eigsMethod_stage2 = GD_plusK;
               }
               else if (strcmp(stringValue, "GD_Olsen_plusK") == 0) {
		       primme_svds->eigsMethod_stage2 = GD_Olsen_plusK;
               }
               else if (strcmp(stringValue, "JD_Olsen_plusK") == 0) {
		       primme_svds->eigsMethod_stage2 = JD_Olsen_plusK;
               }
               else if (strcmp(stringValue, "RQI") == 0) {
		       primme_svds->eigsMethod_stage2 = RQI;
               }
               else if (strcmp(stringValue, "JDQR") == 0) {
		       primme_svds->eigsMethod_stage2 = JDQR;
               }
               else if (strcmp(stringValue, "JDQMR") == 0) {
		       primme_svds->eigsMethod_stage2 = JDQMR;
               }
               else if (strcmp(stringValue, "JDQMR_ETol") == 0) {
		       primme_svds->eigsMethod_stage2 = JDQMR_ETol;
               }
               else if (strcmp(stringValue, "SUBSPACE_ITERATION") == 0) {
		       primme_svds->eigsMethod_stage2 = SUBSPACE_ITERATION;
               }
               else if (strcmp(stringValue, "LOBPCG_OrthoBasis") == 0) {
		       primme_svds->eigsMethod_stage2 = LOBPCG_OrthoBasis;
               }
               else if (strcmp(stringValue, "LOBPCG_OrthoBasis_Window") == 0) {
		       primme_svds->eigsMethod_stage2 = LOBPCG_OrthoBasis_Window;
               }
               else {
                  printf("Invalid eigenMethod value\n");
                  ret = 0;
               }
            }
         }
         else if (strcmp(ident, "primme_svds.numSvals") == 0) {
            ret = fscanf(configFile, "%d\n", &primme_svds->numSvals);
         }
         else if (strcmp(ident, "primme_svds.numSvals") == 0) {
            ret = fscanf(configFile, "%d\n", &primme_svds->numSvals);
         }
	     else if (strcmp(ident, "primme_svds.eps") == 0) {
            ret = fscanf(configFile, "%le\n", &primme_svds->eps);
         }
         else if (strcmp(ident, "primme_svds.aNorm") == 0) {
            ret = fscanf(configFile, "%le\n", &primme_svds->aNorm);
         }
         else if (strcmp(ident, "primme_svds.maxBasisSize") == 0) {
            ret = fscanf(configFile, "%d\n", &primme_svds->maxBasisSize);
         }
         else if (strcmp(ident, "primme_svds.maxBlockSize") == 0) {
            ret = fscanf(configFile, "%d\n", &primme_svds->maxBlockSize);
         }
         else if (strcmp(ident, "primme_svds.initSize") == 0) {
            ret = fscanf(configFile, "%d\n", &primme_svds->initSize);
         }
         else if (strcmp(ident, "primme_svds.maxMatvecs") == 0) {
            ret = fscanf(configFile, "%d\n", &primme_svds->maxMatvecs);
         }
         else if (strcmp(ident, "primme_svds.printLevel") == 0) {
            ret = fscanf(configFile, "%d\n", &primme_svds->printLevel);
         }
         else if (strcmp(ident, "primme_svds.target") == 0) {
            ret = fscanf(configFile, "%s\n", stringValue);
            if (ret == 1) {
               if (strcmp(stringValue, "primme_svds_smallest") == 0) {
                  primme_svds->target = primme_svds_smallest;
               }
               else if (strcmp(stringValue, "primme_svds_largest") == 0) {
                  primme_svds->target = primme_svds_largest;
               }
               else {
                  printf("Invalid target value\n");
                  ret = 0;
               }
            }
         }
         else if (strcmp(ident, "primme_svds.svdsMethod") == 0) {
            ret = fscanf(configFile, "%s\n", stringValue);
            if (ret == 1) {
               if (strcmp(stringValue, "primme_svds_hybrid") == 0) {
                  primme_svds->svdsMethod = primme_svds_hybrid;
               }
               else if (strcmp(stringValue, "primme_svds_normalequations") == 0) {
                  primme_svds->svdsMethod = primme_svds_normalequations;
               }
               else if (strcmp(stringValue, "primme_svds_augmented") == 0) {
                  primme_svds->svdsMethod = primme_svds_augmented;
               }
               else {
                  printf("Invalid svdsMethod value\n");
                  ret = 0;
               }
            }
         }
         else if (strcmp(ident, "primme_svds.precondition") == 0) {
            ret = fscanf(configFile, "%d\n", &primme_svds->precondition);
         }
         else if (strcmp(ident, "primme_svds.iseed") == 0) {                                                                         
            ret = 1;
            for (i=0;i<4; i++) {
                ret = fscanf(configFile, "%d", &primme_svds->iseed[i]);
                if (ret != 1) break;
            }
            if (ret == 1) {
                fgets(ident, 2048, configFile);
            }
         }
         else if (strcmp(ident, "primme_svds.numTargetShifts") == 0) {
            ret = fscanf(configFile, "%d\n", &primme_svds->numTargetShifts);
         }
         else if (strcmp(ident, "primme_svds.targetShifts") == 0) {  
            ret = 1;
            if (primme_svds->numTargetShifts >0) {
                primme_svds->targetShifts = (double *)primme_calloc(
                    primme_svds->numTargetShifts, sizeof(double), "targetShifts");
                for (i=0;i<primme_svds->numTargetShifts; i++) {
                    ret = fscanf(configFile, "%le", &primme_svds->targetShifts[i]);
                    if (ret != 1) break;
                }
            }
            if (ret == 1) {
                fgets(ident, 2048, configFile);
            }
         }
         else {
            fprintf(stderr, 
	       "ERROR(read_solver_params): Invalid parameter '%s'\n", ident);
            return(-1);
         }

         line++;
      }
      else {
         fprintf(stderr, 
	    "ERROR(read_solver_params): Invalid operator on %d\n", line);
         return(-1);
      }

      if (ret != 1) {
         fprintf(stderr, 
	 "ERROR(read_solver_params): Could not read value on line %d\n", line);
         return(-1);
      }
   }

   // Set up the output file in primme_svds, from the filename read in driverConfig
   if (primme_svds->procID == 0) {
      if (strcmp(outputFileName, "stdout") != 0) {
         if ((primme_svds->outputFile = fopen(outputFileName, "w+")) == NULL) {
            fprintf(stderr, 
		   "ERROR(read_solver_params): Could not open output file\n");
         }
      }
      else {
         primme_svds->outputFile = stdout;
      }
   }
   else {
      primme_svds->outputFile = stdout;
   }

   return (0);
}

/******************************************************************************
 *
 * Reads the parameters necessary for the test driver
 * eg., matrix, preconditioning parameters and choice, output files, etc
 * This function does not read any solver parameters.
 *
******************************************************************************/
int read_driver_params(char *configFileName, driver_params *driver) {

   int line, ret;
   char ident[128];
   char op[128];
   FILE *configFile;


   if ((configFile = fopen(configFileName, "r")) == NULL) {
      fprintf(stderr,"Error(read_driver_params): Could not open config file\n");
      fprintf(stderr,"Driver config file: %s\n", configFileName);
      return(-1);
   }

   line = 1;
   while (EOF != fscanf(configFile, "%s", ident)) {
      if (strcmp(ident, "//") == 0) {
         fgets(ident, 2048, configFile);
         line++;
         continue;
      }
      else {
         ret = fscanf(configFile, "%s", op);
      }

      if (ret == EOF) {
         fprintf(stderr, "ERROR(read_driver_params): Unexpected end of file\n");
         return(-1);
      }

      if (strcmp(op, "=") == 0) {
	 // Matrix, partitioning and I/O params 
         if (strcmp(ident, "driver.outputFile") == 0) {
            ret = fscanf(configFile, "%s\n", driver->outputFileName);
         }
         else if (strcmp(ident, "driver.partId") == 0) {
            ret = fscanf(configFile, "%s\n", driver->partId);
         }
         else if (strcmp(ident, "driver.partDir") == 0) {
            ret = fscanf(configFile, "%s\n", driver->partDir);
         }
         else if (strcmp(ident, "driver.matrixFile") == 0) {
            ret = fscanf(configFile, "%s\n", driver->matrixFileName);
         }
	 // Preconditioning parameters
         else if (strcmp(ident, "driver.PrecChoice") == 0) {
            ret = fscanf(configFile, "%d\n", &driver->PrecChoice);
         }
         else if (strcmp(ident, "driver.shift") == 0) {
            ret = fscanf(configFile, "%le\n", &driver->shift);
         }
         else if (strcmp(ident, "driver.isymm") == 0) {
            ret = fscanf(configFile, "%d\n", &driver->isymm);
         }
         else if (strcmp(ident, "driver.level") == 0) {
            ret = fscanf(configFile, "%d\n", &driver->level);
         }
         else if (strcmp(ident, "driver.threshold") == 0) {
            ret = fscanf(configFile, "%lf\n", &driver->threshold);
         }
         else if (strcmp(ident, "driver.filter") == 0) {
            ret = fscanf(configFile, "%lf\n", &driver->filter);
	 }
         else {
            fprintf(stderr, 
	      "ERROR(read_driver_params): Invalid parameter '%s'\n", ident);
            return(-1);
         }

         line++;
      }
      else {
         fprintf(stderr, "ERROR(read_driver_params): Invalid operator on %d\n",
		 line);
         return(-1);
      }

      if (ret != 1) {
         fprintf(stderr, 
	  "ERROR(read_driver_params): Could not read value on line %d\n", line);
         return(-1);
      }
   }
   fclose(configFile);
   return (0);
}


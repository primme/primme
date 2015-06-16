/*  PRIMME PReconditioned Iterative MultiMethod Eigensolver
 *   Copyright (C) 2005  James R. McCombs,  Andreas Stathopoulos
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
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include "shared_utils.h"


/******************************************************************************
 *
 * Reads the solver parameters for configuring the primme structure of 
 * PRIMME. "method" is not part of the struct primme, but if specified
 * primme will be setup accordingly by primme_set_method().
 *
******************************************************************************/
int read_solver_params(char *configFileName, char *outputFileName,
		primme_params *primme, primme_preset_method *method) {

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

   *method = (primme_preset_method) -1;   // Set it as undefined

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
         if (strcmp(ident, "method") == 0) {
            ret = fscanf(configFile, "%s\n", stringValue);
            if (ret == 1) {
               if (strcmp(stringValue,      "DYNAMIC") == 0) {
		       *method = DYNAMIC;
               }
	       else if (strcmp(stringValue, "DEFAULT_MIN_TIME") == 0) {
		       *method = DEFAULT_MIN_TIME;
               }
               else if (strcmp(stringValue, "DEFAULT_MIN_MATVECS") == 0) {
		       *method = DEFAULT_MIN_MATVECS;
               }
               else if (strcmp(stringValue, "Arnoldi") == 0) {
		       *method = Arnoldi;
               }
               else if (strcmp(stringValue, "GD") == 0) {
		       *method = GD;
               }
               else if (strcmp(stringValue, "GD_plusK") == 0) {
		       *method = GD_plusK;
               }
               else if (strcmp(stringValue, "GD_Olsen_plusK") == 0) {
		       *method = GD_Olsen_plusK;
               }
               else if (strcmp(stringValue, "JD_Olsen_plusK") == 0) {
		       *method = JD_Olsen_plusK;
               }
               else if (strcmp(stringValue, "RQI") == 0) {
		       *method = RQI;
               }
               else if (strcmp(stringValue, "JDQR") == 0) {
		       *method = JDQR;
               }
               else if (strcmp(stringValue, "JDQMR") == 0) {
		       *method = JDQMR;
               }
               else if (strcmp(stringValue, "JDQMR_ETol") == 0) {
		       *method = JDQMR_ETol;
               }
               else if (strcmp(stringValue, "SUBSPACE_ITERATION") == 0) {
		       *method = SUBSPACE_ITERATION;
               }
               else if (strcmp(stringValue, "LOBPCG_OrthoBasis") == 0) {
		       *method = LOBPCG_OrthoBasis;
               }
               else if (strcmp(stringValue, "LOBPCG_OrthoBasis_Window") == 0) {
		       *method = LOBPCG_OrthoBasis_Window;
               }
               else {
                  printf("Invalid target value\n");
                  ret = 0;
               }
            }
         }
         else if (strcmp(ident, "primme.numEvals") == 0) {
            ret = fscanf(configFile, "%d\n", &primme->numEvals);
         }
	 else if (strcmp(ident, "primme.eps") == 0) {
            ret = fscanf(configFile, "%le\n", &primme->eps);
         }
         else if (strcmp(ident, "primme.aNorm") == 0) {
            ret = fscanf(configFile, "%le\n", &primme->aNorm);
         }
         else if (strcmp(ident, "primme.maxBasisSize") == 0) {
            ret = fscanf(configFile, "%d\n", &primme->maxBasisSize);
         }
         else if (strcmp(ident, "primme.minRestartSize") == 0) {
            ret = fscanf(configFile, "%d\n", &primme->minRestartSize);
         }
         else if (strcmp(ident, "primme.locking") == 0) {
            ret = fscanf(configFile, "%d\n", &primme->locking);
         }
         else if (strcmp(ident, "primme.dynamicMethodSwitch") == 0) {
            ret = fscanf(configFile, "%d\n", &primme->dynamicMethodSwitch);
         }
         else if (strcmp(ident, "primme.maxBlockSize") == 0) {
            ret = fscanf(configFile, "%d\n", &primme->maxBlockSize);
         }
         else if (strcmp(ident, "primme.initSize") == 0) {
            ret = fscanf(configFile, "%d\n", &primme->initSize);
         }
         else if (strcmp(ident, "primme.numOrthoConst") == 0) {
            ret = fscanf(configFile, "%d\n", &primme->numOrthoConst);
         }
         else if (strcmp(ident, "primme.maxOuterIterations") == 0) {
            ret = fscanf(configFile, "%d\n", &primme->maxOuterIterations);
         }
         else if (strcmp(ident, "primme.maxMatvecs") == 0) {
            ret = fscanf(configFile, "%d\n", &primme->maxMatvecs);
         }
         else if (strcmp(ident, "primme.printLevel") == 0) {
            ret = fscanf(configFile, "%d\n", &primme->printLevel);
         }
         else if (strcmp(ident, "primme.restarting.scheme") == 0) {
            ret = fscanf(configFile, "%s\n", stringValue); 
            if (ret == 1) {
               if (strcmp(stringValue, "primme_thick") == 0) {
                  primme->restartingParams.scheme = primme_thick;
               }
               else if (strcmp(stringValue, "primme_dtr") == 0) {
                  primme->restartingParams.scheme = primme_dtr;
               }
               else {
                  printf("Invalid restart.scheme value\n");
                  ret = 0;
               }
            }
         }
         else if (strcmp(ident, "primme.restarting.maxPrevRetain") == 0) {
            ret = fscanf(configFile, "%d\n", 
                     &primme->restartingParams.maxPrevRetain);
         }
         else if (strcmp(ident, "primme.target") == 0) {
            ret = fscanf(configFile, "%s\n", stringValue);
            if (ret == 1) {
               if (strcmp(stringValue, "primme_smallest") == 0) {
                  primme->target = primme_smallest;
               }
               else if (strcmp(stringValue, "primme_largest") == 0) {
                  primme->target = primme_largest;
               }
               else if (strcmp(stringValue, "primme_closest_geq") == 0) {
                  primme->target = primme_closest_geq;
               }
               else if (strcmp(stringValue, "primme_closest_leq") == 0) {
                  primme->target = primme_closest_leq;
               }
               else if (strcmp(stringValue, "primme_closest_abs") == 0) {
                  primme->target = primme_closest_abs;
               }
               else {
                  printf("Invalid target value\n");
                  ret = 0;
               }
            }
         }
         else if (strcmp(ident, "primme.numTargetShifts") == 0) {
            ret = fscanf(configFile, "%d\n", &primme->numTargetShifts);
         }
         else if (strcmp(ident, "primme.targetShifts") == 0) {
            ret = 1;
            if (primme->numTargetShifts >0) {
	       primme->targetShifts = (double *)primme_calloc(
	          primme->numTargetShifts, sizeof(double), "targetShifts");
	       for (i=0;i<primme->numTargetShifts; i++) {
                  ret = fscanf(configFile, "%le", &primme->targetShifts[i]);
	          if (ret != 1) break;
	       }
	    }
	    if (ret == 1) {
               fgets(ident, 2048, configFile);
	    }
	 }
         else if (strcmp(ident, "primme.correction.projectors.LeftQ") == 0) {
            ret = fscanf(configFile, "%d\n", 
                     &primme->correctionParams.projectors.LeftQ );
         }
         else if (strcmp(ident, "primme.correction.projectors.LeftX") == 0) {
            ret = fscanf(configFile, "%d\n", 
                     &primme->correctionParams.projectors.LeftX );
         }
         else if (strcmp(ident, "primme.correction.projectors.RightQ") == 0) {
            ret = fscanf(configFile, "%d\n", 
                     &primme->correctionParams.projectors.RightQ );
         }
         else if (strcmp(ident, "primme.correction.projectors.RightX") == 0) {
            ret = fscanf(configFile, "%d\n", 
                     &primme->correctionParams.projectors.RightX );
         }
         else if (strcmp(ident, "primme.correction.projectors.SkewQ") == 0) {
            ret = fscanf(configFile, "%d\n", 
                     &primme->correctionParams.projectors.SkewQ );
         }
         else if (strcmp(ident, "primme.correction.projectors.SkewX") == 0) {
            ret = fscanf(configFile, "%d\n", 
                     &primme->correctionParams.projectors.SkewX );
         }
         else if (strcmp(ident, "primme.correction.convTest") == 0) {
            ret = fscanf(configFile, "%s\n", stringValue);
            if (ret == 1) {
               if (strcmp(stringValue, "primme_full_LTolerance") == 0) {
                  primme->correctionParams.convTest = primme_full_LTolerance;
               }
               else if (strcmp(stringValue, "primme_adaptive") == 0) {
                  primme->correctionParams.convTest = primme_adaptive;
               }
               else if (strcmp(stringValue,"primme_decreasing_LTolerance") == 0)
	       {  primme->correctionParams.convTest = 
		       				 primme_decreasing_LTolerance;
               }
               else if (strcmp(stringValue,"primme_adaptive_ETolerance") == 0) {
                 primme->correctionParams.convTest = primme_adaptive_ETolerance;
              }
            }
         }
         else if (strcmp(ident, "primme.correction.precondition") == 0) {
            ret = fscanf(configFile, "%d\n", 
                         &primme->correctionParams.precondition);
         }
         else if (strcmp(ident, "primme.correction.robustShifts") == 0) {
            ret = fscanf(configFile, "%d\n", 
               &primme->correctionParams.robustShifts);
         }
         else if (strcmp(ident, "primme.correction.maxInnerIterations") == 0) {
            ret = fscanf(configFile, "%d\n", 
               &primme->correctionParams.maxInnerIterations);
         }
         else if (strcmp(ident, "primme.correction.relTolBase") == 0) {
            ret = fscanf(configFile, "%lf\n", 
               &primme->correctionParams.relTolBase);
         }
         else if (strcmp(ident, "primme.iseed") == 0) {
            ret = 1;
	    for (i=0;i<4; i++) {
               ret = fscanf(configFile, "%d", &primme->iseed[i]);
	       if (ret != 1) break;
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

   // Set up the output file in primme, from the filename read in driverConfig
   if (primme->procID == 0) {
      if (strcmp(outputFileName, "stdout") != 0) {
         if ((primme->outputFile = fopen(outputFileName, "w+")) == NULL) {
            fprintf(stderr, 
		   "ERROR(read_solver_params): Could not open output file\n");
         }
      }
      else {
         primme->outputFile = stdout;
      }
   }
   else {
      primme->outputFile = stdout;
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

   return (0);
}


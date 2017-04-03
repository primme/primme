/*******************************************************************************
 *   PRIMME PReconditioned Iterative MultiMethod Eigensolver
 *   Copyright (C) 2017 College of William & Mary,
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
 * File: shared_utils.c
 * 
 * Purpose - Functions to read and print driver_params and primme_params.
 * 
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <assert.h>
#include "shared_utils.h"


/******************************************************************************
 *
 * Reads the solver parameters for configuring the primme structure of 
 * dprimme(). "method" is not part of the struct primme, but if specified
 * primme will be setup accordingly by primme_set_method().
 *
******************************************************************************/
int read_solver_params(char *configFileName, char *outputFileName,
                primme_params *primme, const char* primmeprefix,
                primme_preset_method *method, const char* methodstr) {

   int line, ret, i;
   char ident[2048], *field;
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
      if (strncmp(ident, "//", 2) == 0) {
         if (fgets(ident, 2048, configFile) == NULL) {
            break;
         }
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
         if (strcmp(ident, methodstr) == 0) {
            ret = fscanf(configFile, "%s", stringValue);
            if (ret == 1) {
               ret = 0;
               #define READ_METHOD(V) if (strcmp(stringValue, #V) == 0) {*method = V; ret=1;}
               READ_METHOD(PRIMME_DEFAULT_METHOD);
               READ_METHOD(PRIMME_DYNAMIC);
               READ_METHOD(PRIMME_DEFAULT_MIN_TIME);
               READ_METHOD(PRIMME_DEFAULT_MIN_MATVECS);
               READ_METHOD(PRIMME_Arnoldi);
               READ_METHOD(PRIMME_GD);
               READ_METHOD(PRIMME_GD_plusK);
               READ_METHOD(PRIMME_GD_Olsen_plusK);
               READ_METHOD(PRIMME_JD_Olsen_plusK);
               READ_METHOD(PRIMME_RQI);
               READ_METHOD(PRIMME_JDQR);
               READ_METHOD(PRIMME_JDQMR);
               READ_METHOD(PRIMME_JDQMR_ETol);
               READ_METHOD(PRIMME_STEEPEST_DESCENT);
               READ_METHOD(PRIMME_LOBPCG_OrthoBasis);
               READ_METHOD(PRIMME_LOBPCG_OrthoBasis_Window);
               #undef READ_METHOD
            }
            if (ret == 0) {
               printf("Invalid %s value\n", methodstr);
               return -1;
            }
            line++;
            continue;
         }
         else if (strncmp(ident, primmeprefix, strlen(primmeprefix)) != 0) {
            if (fgets(ident, 2048, configFile) == NULL) {
               break;
            }
            line++;
            continue;
         }
         field = ident + strlen(primmeprefix);

         #define READ_FIELD(V, P) if (strcmp(field, #V) == 0) \
            ret = fscanf(configFile, P, &primme-> V);
         #define READ_FIELD_OP(V, P) if (strcmp(field, #V) == 0) { \
            ret = fscanf(configFile, "%s", stringValue); \
            if (ret == 1) { ret=0; P } \
            if (ret == 0) printf("Invalid " #V " value\n"); \
         }
         #define OPTION(F, V) if (strcmp(stringValue, #V) == 0) { primme-> F = V; ret = 1; }
         #define READ_FIELDParams(S, V, P) if (strcmp(field, #S "." #V) == 0) \
            ret = fscanf(configFile, P, &primme-> S ## Params . V);
         #define READ_FIELD_OPParams(S, V, P) if (strcmp(field, #S "." #V) == 0) { \
            ret = fscanf(configFile, "%s", stringValue); \
            if (ret == 1) { ret=0; P } \
            if (ret == 0) printf("Invalid " #S "." #V " value\n"); \
         }
         #define OPTIONParams(S, F, V) if (strcmp(stringValue, #V) == 0) { primme-> S ## Params . F = V; ret = 1; }
  
         READ_FIELD(printLevel, "%d");
         READ_FIELD(numEvals, "%d");
         READ_FIELD(aNorm, "%le");
         READ_FIELD(eps, "%le");
         READ_FIELD(maxBasisSize, "%d");
         READ_FIELD(minRestartSize, "%d");
         READ_FIELD(maxBlockSize, "%d");
         READ_FIELD(maxOuterIterations, "%" PRIMME_INT_P);
         READ_FIELD(maxMatvecs, "%" PRIMME_INT_P);
         READ_FIELD_OP(target,
            OPTION(target, primme_smallest)
            OPTION(target, primme_largest)
            OPTION(target, primme_closest_geq)
            OPTION(target, primme_closest_leq)
            OPTION(target, primme_closest_abs)
         );
         READ_FIELD_OPParams(projection, projection,
            OPTIONParams(projection, projection, primme_proj_default)
            OPTIONParams(projection, projection, primme_proj_RR)
            OPTIONParams(projection, projection, primme_proj_refined)
            OPTIONParams(projection, projection, primme_proj_harmonic)
         );

         READ_FIELD_OP(initBasisMode,
            OPTION(initBasisMode, primme_init_default)
            OPTION(initBasisMode, primme_init_krylov)
            OPTION(initBasisMode, primme_init_random)
            OPTION(initBasisMode, primme_init_user)
         );

         READ_FIELD(numTargetShifts, "%d");
         if (strcmp(field, "targetShifts") == 0) {
            ret = 1;
            if (primme->numTargetShifts > 0) {
               primme->targetShifts = (double *)primme_calloc(
                  primme->numTargetShifts, sizeof(double), "targetShifts");
               for (i=0; i<primme->numTargetShifts; i++) {
                  ret = fscanf(configFile, "%le", &primme->targetShifts[i]);
                  if (ret != 1) break;
               }
            }
            if (ret == 1) {
               if (fgets(ident, 2048, configFile) == NULL) {
                  break;
               }
            }
         }
 
         READ_FIELD(dynamicMethodSwitch, "%d");
         READ_FIELD(locking, "%d");
         READ_FIELD(initSize, "%d");
         READ_FIELD(numOrthoConst, "%d");

         if (strcmp(field, "iseed") == 0) {
            ret = 1;
            for (i=0;i<4; i++) {
               ret = fscanf(configFile, "%" PRIMME_INT_P, &primme->iseed[i]);
               if (ret != 1) break;
            }
            if (ret == 1) {
               if (fgets(ident, 2048, configFile) == NULL) {
                  break;
               }
            }
         }

         READ_FIELD_OPParams(restarting, scheme,
            OPTIONParams(restarting, scheme, primme_thick)
            OPTIONParams(restarting, scheme, primme_dtr)
         );

         READ_FIELDParams(restarting, maxPrevRetain, "%d");

         READ_FIELDParams(correction, precondition, "%d");
         READ_FIELDParams(correction, robustShifts, "%d");
         READ_FIELDParams(correction, maxInnerIterations, "%d");
         READ_FIELDParams(correction, relTolBase, "%lf");

         READ_FIELD_OPParams(correction, convTest,
            OPTIONParams(correction, convTest, primme_full_LTolerance)
            OPTIONParams(correction, convTest, primme_decreasing_LTolerance)
            OPTIONParams(correction, convTest, primme_adaptive_ETolerance)
            OPTIONParams(correction, convTest, primme_adaptive)
         );

         READ_FIELDParams(correction, projectors.LeftQ , "%d");
         READ_FIELDParams(correction, projectors.LeftX , "%d");
         READ_FIELDParams(correction, projectors.RightQ, "%d");
         READ_FIELDParams(correction, projectors.SkewQ , "%d");
         READ_FIELDParams(correction, projectors.RightX, "%d");
         READ_FIELDParams(correction, projectors.SkewX , "%d");

         if (ret == 0) {
            fprintf(stderr, 
               "ERROR(read_solver_params): Invalid parameter '%s'\n", ident);
            return(-1);
         }
         line++;

         #undef READ_FIELD
         #undef READ_FIELD_OP
         #undef OPTION
         #undef READ_FIELDParams
         #undef READ_FIELD_OPParams
         #undef OPTIONParams
      }
      else {
         fprintf(stderr, 
            "ERROR(read_solver_params): Invalid operator on line %d\n", line);
         return(-1);
      }

      if (ret != 1) {
         fprintf(stderr, 
         "ERROR(read_solver_params): Could not read value on line %d\n", line);
         return(-1);
      }
   }

   /* Set up the output file in primme, from the filename read in driverConfig */
   if (primme->procID == 0) {
      if (outputFileName[0] && strcmp(outputFileName, "stdout") != 0) {
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

   fclose(configFile);
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
   char ident[2048];
   char op[128];
   char stringValue[128];
   FILE *configFile;


   memset(driver, 0, sizeof(*driver));
   if ((configFile = fopen(configFileName, "r")) == NULL) {
      fprintf(stderr,"Error(read_driver_params): Could not open config file\n");
      fprintf(stderr,"Driver config file: %s\n", configFileName);
      return(-1);
   }

   line = 1;
   while (EOF != fscanf(configFile, "%s", ident)) {
      if (strncmp(ident, "//", 2) == 0) {
         if (fgets(ident, 2048, configFile) == NULL) {
            break;
         }
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
         /* Matrix, partitioning and I/O params  */
         if (strcmp(ident, "driver.outputFile") == 0) {
            ret = fscanf(configFile, "%s", driver->outputFileName);
         }
         else if (strcmp(ident, "driver.partId") == 0) {
            ret = fscanf(configFile, "%s", driver->partId);
         }
         else if (strcmp(ident, "driver.partDir") == 0) {
            ret = fscanf(configFile, "%s", driver->partDir);
         }
         else if (strcmp(ident, "driver.matrixFile") == 0) {
            ret = fscanf(configFile, "%s", driver->matrixFileName);
         }
         else if (strcmp(ident, "driver.initialGuessesFile") == 0) {
            ret = fscanf(configFile, "%s", driver->initialGuessesFileName);
         }
         else if (strcmp(ident, "driver.initialGuessesPert") == 0) {
            ret = fscanf(configFile, "%le", &driver->initialGuessesPert);
         }
         else if (strcmp(ident, "driver.saveXFile") == 0) {
            ret = fscanf(configFile, "%s", driver->saveXFileName);
         }
         else if (strcmp(ident, "driver.checkXFile") == 0) {
            ret = fscanf(configFile, "%s", driver->checkXFileName);
         }
         else if (strcmp(ident, "driver.checkInterface") == 0) {
            ret = fscanf(configFile, "%d", &driver->checkInterface);
         }
         else if (strcmp(ident, "driver.matrixChoice") == 0) {
            ret = fscanf(configFile, "%s", stringValue);
            if (ret == 1) {
               if (strcmp(stringValue, "default") == 0) {
                  driver->matrixChoice = driver_default;
               }
               else if (strcmp(stringValue, "native") == 0) {
                  driver->matrixChoice = driver_native;
               }
               else if (strcmp(stringValue, "parasails") == 0) {
                  driver->matrixChoice = driver_parasails;
               }
               else if (strcmp(stringValue, "petsc") == 0) {
                  driver->matrixChoice = driver_petsc;
               }
               else if (strcmp(stringValue, "rsb") == 0) {
                  driver->matrixChoice = driver_rsb;
               }
               else {
                  fprintf(stderr, 
                     "ERROR(read_driver_params): Invalid parameter '%s'\n", ident);
               }
            }
         }
         /* Preconditioning parameters */
         else if (strcmp(ident, "driver.PrecChoice") == 0) {
            ret = fscanf(configFile, "%s", stringValue);
            if (ret == 1) {
               if (strcmp(stringValue, "noprecond") == 0) {
                  driver->PrecChoice = driver_noprecond;
               }
               else if (strcmp(stringValue, "jacobi") == 0) {
                  driver->PrecChoice = driver_jacobi;
               }
               else if (strcmp(stringValue, "davidsonjacobi") == 0) {
                  driver->PrecChoice = driver_jacobi_i;
               }
               else if (strcmp(stringValue, "ilut") == 0) {
                  driver->PrecChoice = driver_ilut;
               }
               else if (strcmp(stringValue, "normal") == 0) {
                  driver->PrecChoice = driver_normal;
               }
               else if (strcmp(stringValue, "bjacobi") == 0) {
                  driver->PrecChoice = driver_bjacobi;
               }
               else {
                  fprintf(stderr, 
                     "ERROR(read_driver_params): Invalid parameter '%s'\n", ident);
               }
            }
         }
         else if (strcmp(ident, "driver.shift") == 0) {
            ret = fscanf(configFile, "%le", &driver->shift);
         }
         else if (strcmp(ident, "driver.isymm") == 0) {
            ret = fscanf(configFile, "%d", &driver->isymm);
         }
         else if (strcmp(ident, "driver.level") == 0) {
            ret = fscanf(configFile, "%d", &driver->level);
         }
         else if (strcmp(ident, "driver.threshold") == 0) {
            ret = fscanf(configFile, "%lf", &driver->threshold);
         }
         else if (strcmp(ident, "driver.filter") == 0) {
            ret = fscanf(configFile, "%lf", &driver->filter);
         }
         else if (strncmp(ident, "driver.", 7) == 0) {
            fprintf(stderr, 
              "ERROR(read_driver_params): Invalid parameter '%s'\n", ident);
            return(-1);
         }
         else {
            if (fgets(ident, 2048, configFile) == NULL) {
               break;
            }
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

void driver_display_params(driver_params driver, FILE *outputFile) {

const char *strPrecChoice[] = {"noprecond", "jacobi", "davidsonjacobi", "ilut", "normal", "bjacobi"};
const char *strMatrixChoice[] = {"default", "native", "petsc", "parasails", "rsb"};
 
fprintf(outputFile, "// ---------------------------------------------------\n"
                    "//                 driver configuration               \n"
                    "// ---------------------------------------------------\n");

fprintf(outputFile, "driver.partId        = %s\n", driver.partId);
fprintf(outputFile, "driver.partDir       = %s\n", driver.partDir);
fprintf(outputFile, "driver.matrixFile    = %s\n", driver.matrixFileName);
fprintf(outputFile, "driver.matrixChoice  = %s\n", strMatrixChoice[driver.matrixChoice]);
fprintf(outputFile, "driver.initialGuessesFile = %s\n", driver.initialGuessesFileName);
fprintf(outputFile, "driver.initialGuessesPert = %e\n", driver.initialGuessesPert);
fprintf(outputFile, "driver.saveXFile     = %s\n", driver.saveXFileName);
fprintf(outputFile, "driver.checkXFile    = %s\n", driver.checkXFileName);
fprintf(outputFile, "driver.checkInterface = %d\n", driver.checkInterface);
fprintf(outputFile, "driver.PrecChoice    = %s\n", strPrecChoice[driver.PrecChoice]);
fprintf(outputFile, "driver.shift         = %e\n", driver.shift);
fprintf(outputFile, "driver.isymm         = %d\n", driver.isymm);
fprintf(outputFile, "driver.level         = %d\n", driver.level);
fprintf(outputFile, "driver.threshold     = %f\n", driver.threshold);
fprintf(outputFile, "driver.filter        = %f\n\n", driver.filter);

}

void driver_display_method(primme_preset_method method, const char* methodstr, FILE *outputFile) {

   const char *strMethod[] = {
      "PRIMME_DEFAULT_METHOD",
      "PRIMME_DYNAMIC",
      "PRIMME_DEFAULT_MIN_TIME",
      "PRIMME_DEFAULT_MIN_MATVECS",
      "PRIMME_Arnoldi",
      "PRIMME_GD",
      "PRIMME_GD_plusK",
      "PRIMME_GD_Olsen_plusK",
      "PRIMME_JD_Olsen_plusK",
      "PRIMME_RQI",
      "PRIMME_JDQR",
      "PRIMME_JDQMR",
      "PRIMME_JDQMR_ETol",
      "PRIMME_STEEPEST_DESCENT",
      "PRIMME_LOBPCG_OrthoBasis",
      "PRIMME_LOBPCG_OrthoBasis_Window"};

   fprintf(outputFile, "%s               = %s\n", methodstr, strMethod[method]);

}

void driver_display_methodsvd(primme_svds_preset_method method, const char* methodstr, FILE *outputFile) {

   const char *strMethod[] = {
      "primme_svds_default",
      "primme_svds_hybrid",
      "primme_svds_normalequations",
      "primme_svds_augmented"};

   fprintf(outputFile, "%s               = %s\n", methodstr, strMethod[method]);

}


int read_solver_params_svds(char *configFileName, char *outputFileName,
                primme_svds_params *primme_svds, const char* primmeprefix,
                primme_svds_preset_method *method, const char* methodstr,
                primme_preset_method *primme_method,
                primme_preset_method *primme_methodStage2) {

   int line, ret, i;
   char ident[2048], *field;
   char op[128];
   char stringValue[128];
   FILE *configFile;

   if ((configFile = fopen(configFileName, "r")) == NULL) {
      fprintf(stderr,"Error(read_solver_params_svds): Could not open config file\n");
      fprintf(stderr, "config file: %s\n", configFileName);
      return(-1);
   }

   line = 1;
   while (EOF != fscanf(configFile, "%s", ident)) {
      if (strncmp(ident, "//", 2) == 0) {
         if (fgets(ident, 2048, configFile) == NULL) {
            break;
         }
         line++;
         continue;
      }
      else {
         ret = fscanf(configFile, "%s", op);
      }

      if (ret == EOF) {
         fprintf(stderr, "ERROR(read_solver_params_svds): Unexpected end of file\n");
         return(-1);
      }

      if (strcmp(op, "=") == 0) {
         if (strcmp(ident, methodstr) == 0) {
            ret = fscanf(configFile, "%s", stringValue);
            if (ret == 1) {
               ret = 0;
               #define READ_METHOD(V) if (strcmp(stringValue, #V) == 0) {*method = V; ret=1;}
               READ_METHOD(primme_svds_default);
               READ_METHOD(primme_svds_hybrid);
               READ_METHOD(primme_svds_normalequations);
               READ_METHOD(primme_svds_augmented);
               #undef READ_METHOD
            }
            if (ret == 0) {
               printf("Invalid %s value\n", methodstr);
               return -1;
            }
            line++;
            continue;
         }
         else if (strncmp(ident, primmeprefix, strlen(primmeprefix)) != 0) {
            if (fgets(ident, 2048, configFile) == NULL) {
               break;
            }
            line++;
            continue;
         }
         field = ident + strlen(primmeprefix);

         #define READ_FIELD(V, P) if (strcmp(field, #V) == 0) \
            ret = fscanf(configFile, P, &primme_svds-> V);
         #define READ_FIELD_OP(V, P) if (strcmp(field, #V) == 0) { \
            ret = fscanf(configFile, "%s", stringValue); \
            if (ret == 1) { ret=0; P } \
            if (ret == 0) printf("Invalid " #V " value\n"); \
         }
         #define OPTION(F, V) if (strcmp(stringValue, #V) == 0) { primme_svds-> F = V; ret = 1; }
  
         READ_FIELD(printLevel, "%d");
         READ_FIELD(numSvals, "%d");
         READ_FIELD(aNorm, "%le");
         READ_FIELD(eps, "%le");
         READ_FIELD(maxBasisSize, "%d");
         READ_FIELD(maxBlockSize, "%d");
         READ_FIELD(maxMatvecs, "%" PRIMME_INT_P);

         READ_FIELD_OP(target,
            OPTION(target, primme_svds_smallest)
            OPTION(target, primme_svds_largest)
            OPTION(target, primme_svds_closest_abs)
         );

         READ_FIELD(numTargetShifts, "%d");
         if (strcmp(field, "targetShifts") == 0) {
            ret = 1;
            if (primme_svds->numTargetShifts > 0) {
               primme_svds->targetShifts = (double *)primme_calloc(
                  primme_svds->numTargetShifts, sizeof(double), "targetShifts");
               for (i=0; i<primme_svds->numTargetShifts; i++) {
                  ret = fscanf(configFile, "%le", &primme_svds->targetShifts[i]);
                  if (ret != 1) break;
               }
            }
            if (ret == 1) {
               if (fgets(ident, 2048, configFile) == NULL) {
                  break;
               }
            }
         }
 
         READ_FIELD(locking, "%d");
         READ_FIELD(initSize, "%d");
         READ_FIELD(numOrthoConst, "%d");

         if (strcmp(field, "iseed") == 0) {
            ret = 1;
            for (i=0;i<4; i++) {
               ret = fscanf(configFile, "%" PRIMME_INT_P, &primme_svds->iseed[i]);
               if (ret != 1) break;
            }
            if (ret == 1) {
               if (fgets(ident, 2048, configFile) == NULL) {
                  break;
               }
            }
         }

         READ_FIELD(precondition, "%d");

         READ_FIELD_OP(method,
            OPTION(method, primme_svds_op_none)
            OPTION(method, primme_svds_op_AtA)
            OPTION(method, primme_svds_op_AAt)
            OPTION(method, primme_svds_op_augmented)
         );

         READ_FIELD_OP(methodStage2,
            OPTION(method, primme_svds_op_none)
            OPTION(method, primme_svds_op_AtA)
            OPTION(method, primme_svds_op_AAt)
            OPTION(method, primme_svds_op_augmented)
         );

         if (ret == 0) {
            fprintf(stderr, 
               "ERROR(read_solver_params_svds): Invalid parameter '%s'\n", ident);
            return(-1);
         }
         line++;

         #undef READ_FIELD
         #undef READ_FIELD_OP
         #undef OPTION
       }
      else {
         fprintf(stderr, 
            "ERROR(read_solver_params_svds): Invalid operator on %d\n", line);
         return(-1);
      }

      if (ret != 1) {
         fprintf(stderr, 
         "ERROR(read_solver_params_svds): Could not read value on line %d\n", line);
         return(-1);
      }
   }

   /* Set up the output file in primme_svds, from the filename read in driverConfig */
   if (primme_svds->procID == 0) {
      if (outputFileName[0] && strcmp(outputFileName, "stdout") != 0) {
         if ((primme_svds->outputFile = fopen(outputFileName, "w+")) == NULL) {
            fprintf(stderr, 
                   "ERROR(read_solver_params_svds): Could not open output file\n");
         }
      }
      else {
         primme_svds->outputFile = stdout;
      }
   }
   else {
      primme_svds->outputFile = stdout;
   }

   fclose(configFile);

   read_solver_params(configFileName, outputFileName, &primme_svds->primme,
                      "primme.", primme_method, "primme.method");
   read_solver_params(configFileName, outputFileName, &primme_svds->primmeStage2,
                      "primmeStage2.", primme_methodStage2, "primmeStage2.method");

   return (0);
}

#ifdef USE_MPI

/******************************************************************************
 * Function to broadcast the primme data structure to all processors
 *
 * EXCEPTIONS: procID and seed[] are not copied from processor 0. 
 *             Each process creates their own.
******************************************************************************/
void broadCast(primme_params *primme, primme_preset_method *method, 
   driver_params *driver, int master, MPI_Comm comm){

   int i;

   if (driver) {
      MPI_Bcast(driver->outputFileName, 512, MPI_CHAR, 0, comm);
      MPI_Bcast(driver->matrixFileName, 1024, MPI_CHAR, 0, comm);
      MPI_Bcast(driver->initialGuessesFileName, 1024, MPI_CHAR, 0, comm);
      MPI_Bcast(driver->saveXFileName, 1024, MPI_CHAR, 0, comm);
      MPI_Bcast(driver->checkXFileName, 1024, MPI_CHAR, 0, comm);
      MPI_Bcast(&driver->initialGuessesPert, 1, MPI_DOUBLE, 0, comm);
      MPI_Bcast(&driver->matrixChoice, 1, MPI_INT, 0, comm);
      MPI_Bcast(&driver->PrecChoice, 1, MPI_INT, 0, comm);
      MPI_Bcast(&driver->isymm, 1, MPI_INT, 0, comm);
      MPI_Bcast(&driver->level, 1, MPI_INT, 0, comm);
      MPI_Bcast(&driver->threshold, 1, MPI_DOUBLE, 0, comm);
      MPI_Bcast(&driver->filter, 1, MPI_DOUBLE, 0, comm);
      MPI_Bcast(&driver->shift, 1, MPI_DOUBLE, 0, comm);
   }

   MPI_Bcast(&(primme->numEvals), 1, MPI_INT, 0, comm);
   MPI_Bcast(&(primme->target), 1, MPI_INT, 0, comm);
   MPI_Bcast(&(primme->numTargetShifts), 1, MPI_INT, 0, comm);

   if (primme->numTargetShifts > 0 && !master) {
      primme->targetShifts = (double *)primme_calloc(
         primme->numTargetShifts, sizeof(double), "targetShifts");
   }
   assert(!master || primme->numTargetShifts == 0 || primme->targetShifts);
   MPI_Bcast(primme->targetShifts, primme->numTargetShifts, MPI_DOUBLE, 0, comm);

   MPI_Bcast(&(primme->locking), 1, MPI_INT, 0, comm);
   MPI_Bcast(&(primme->dynamicMethodSwitch), 1, MPI_INT, 0, comm);
   MPI_Bcast(&(primme->initSize), 1, MPI_INT, 0, comm);
   MPI_Bcast(&(primme->numOrthoConst), 1, MPI_INT, 0, comm);
   MPI_Bcast(&(primme->maxBasisSize), 1, MPI_INT, 0, comm);
   MPI_Bcast(&(primme->minRestartSize), 1, MPI_INT, 0, comm);
   MPI_Bcast(&(primme->maxBlockSize), 1, MPI_INT, 0, comm);
   MPI_Bcast(&(primme->maxMatvecs), 1, MPI_INT, 0, comm);
   MPI_Bcast(&(primme->maxOuterIterations), 1, MPI_INT, 0, comm);
   MPI_Bcast(&(primme->aNorm), 1, MPI_DOUBLE, 0, comm);
   MPI_Bcast(&(primme->eps), 1, MPI_DOUBLE, 0, comm);
   MPI_Bcast(&(primme->printLevel), 1, MPI_INT, 0, comm);
   MPI_Bcast(&(primme->initBasisMode), 1, MPI_INT, 0, comm);

   MPI_Bcast(&(primme->projectionParams.projection), 1, MPI_INT, 0, comm);
   MPI_Bcast(&(primme->restartingParams.scheme), 1, MPI_INT, 0, comm);
   MPI_Bcast(&(primme->restartingParams.maxPrevRetain), 1, MPI_INT, 0, comm);

   MPI_Bcast(&(primme->correctionParams.precondition), 1, MPI_INT, 0, comm);
   MPI_Bcast(&(primme->correctionParams.robustShifts), 1, MPI_INT, 0, comm);
   MPI_Bcast(&(primme->correctionParams.maxInnerIterations),1, MPI_INT, 0,comm);
   MPI_Bcast(&(primme->correctionParams.convTest), 1, MPI_INT, 0, comm);
   MPI_Bcast(&(primme->correctionParams.relTolBase), 1, MPI_DOUBLE, 0, comm);
   MPI_Bcast(&(primme->correctionParams.projectors.LeftQ),  1, MPI_INT, 0,comm);
   MPI_Bcast(&(primme->correctionParams.projectors.LeftX),  1, MPI_INT, 0,comm);
   MPI_Bcast(&(primme->correctionParams.projectors.RightQ), 1, MPI_INT, 0,comm);
   MPI_Bcast(&(primme->correctionParams.projectors.RightX), 1, MPI_INT, 0,comm);
   MPI_Bcast(&(primme->correctionParams.projectors.SkewQ),  1, MPI_INT, 0,comm);
   MPI_Bcast(&(primme->correctionParams.projectors.SkewX),  1, MPI_INT, 0,comm);
   MPI_Bcast(&(primme->correctionParams.projectors.SkewX),  1, MPI_INT, 0,comm);

   MPI_Bcast(method, 1, MPI_INT, 0, comm);
}

/******************************************************************************
 * Function to broadcast the primme svds data structure to all processors
 *
 * EXCEPTIONS: procID and seed[] are not copied from processor 0. 
 *             Each process creates their own.
******************************************************************************/
void broadCast_svds(primme_svds_params *primme_svds, primme_svds_preset_method *method,
   primme_preset_method *primmemethod, primme_preset_method *primmemethodStage2,
   driver_params *driver, int master, MPI_Comm comm){

   MPI_Bcast(driver->outputFileName, 512, MPI_CHAR, 0, comm);
   MPI_Bcast(driver->matrixFileName, 1024, MPI_CHAR, 0, comm);
   MPI_Bcast(driver->initialGuessesFileName, 1024, MPI_CHAR, 0, comm);
   MPI_Bcast(driver->saveXFileName, 1024, MPI_CHAR, 0, comm);
   MPI_Bcast(driver->checkXFileName, 1024, MPI_CHAR, 0, comm);
   MPI_Bcast(&driver->initialGuessesPert, 1, MPI_DOUBLE, 0, comm);
   MPI_Bcast(&driver->matrixChoice, 1, MPI_INT, 0, comm);
   MPI_Bcast(&driver->PrecChoice, 1, MPI_INT, 0, comm);
   MPI_Bcast(&driver->isymm, 1, MPI_INT, 0, comm);
   MPI_Bcast(&driver->level, 1, MPI_INT, 0, comm);
   MPI_Bcast(&driver->threshold, 1, MPI_DOUBLE, 0, comm);
   MPI_Bcast(&driver->filter, 1, MPI_DOUBLE, 0, comm);
   MPI_Bcast(&driver->shift, 1, MPI_DOUBLE, 0, comm);

   MPI_Bcast(&(primme_svds->numSvals), 1, MPI_INT, 0, comm);
   MPI_Bcast(&(primme_svds->target), 1, MPI_INT, 0, comm);
   MPI_Bcast(&(primme_svds->numTargetShifts), 1, MPI_INT, 0, comm);

   if (primme_svds->numTargetShifts > 0 && !master) {
      primme_svds->targetShifts = (double *)primme_calloc(
         primme_svds->numTargetShifts, sizeof(double), "targetShifts");
   }
   assert(!master || primme_svds->numTargetShifts == 0 || primme_svds->targetShifts);
   MPI_Bcast(primme_svds->targetShifts, primme_svds->numTargetShifts, MPI_DOUBLE, 0, comm);
   MPI_Bcast(&(primme_svds->locking), 1, MPI_INT, 0, comm);
   MPI_Bcast(&(primme_svds->initSize), 1, MPI_INT, 0, comm);
   MPI_Bcast(&(primme_svds->numOrthoConst), 1, MPI_INT, 0, comm);
   MPI_Bcast(&(primme_svds->maxBasisSize), 1, MPI_INT, 0, comm);
   MPI_Bcast(&(primme_svds->maxBlockSize), 1, MPI_INT, 0, comm);
   MPI_Bcast(&(primme_svds->maxMatvecs), 1, MPI_INT, 0, comm);
   MPI_Bcast(&(primme_svds->aNorm), 1, MPI_DOUBLE, 0, comm);
   MPI_Bcast(&(primme_svds->eps), 1, MPI_DOUBLE, 0, comm);
   MPI_Bcast(&(primme_svds->printLevel), 1, MPI_INT, 0, comm);
   MPI_Bcast(&(primme_svds->method), 1, MPI_INT, 0, comm);
   MPI_Bcast(&(primme_svds->methodStage2), 1, MPI_INT, 0, comm);
   MPI_Bcast(&(primme_svds->precondition), 1, MPI_INT, 0, comm);

   MPI_Bcast(method, 1, MPI_INT, 0, comm);
   broadCast(&primme_svds->primme, primmemethod,  NULL, master, comm);
   broadCast(&primme_svds->primmeStage2, primmemethodStage2,  NULL, master, comm);
}

#ifdef USE_PETSC
#include <petscmat.h>
#endif

/******************************************************************************
 * MPI globalSumDouble function
 *
******************************************************************************/
void par_GlobalSumDouble(void *sendBuf, void *recvBuf, int *count, 
                         primme_params *primme, int *ierr) {
   MPI_Comm communicator = *(MPI_Comm *) primme->commInfo;

#ifdef USE_PETSC
   extern PetscLogEvent PRIMME_GLOBAL_SUM;
   PetscLogEventBegin(PRIMME_GLOBAL_SUM,0,0,0,0);
   *ierr = MPI_Allreduce(sendBuf, recvBuf, *count, MPI_DOUBLE, MPI_SUM, communicator);
   PetscLogEventEnd(PRIMME_GLOBAL_SUM,0,0,0,0);
#else
   *ierr = MPI_Allreduce(sendBuf, recvBuf, *count, MPI_DOUBLE, MPI_SUM, communicator);
#endif
}

void par_GlobalSumDoubleSvds(void *sendBuf, void *recvBuf, int *count, 
                         primme_svds_params *primme_svds, int *ierr) {
   MPI_Comm communicator = *(MPI_Comm *) primme_svds->commInfo;

#ifdef USE_PETSC
   extern PetscLogEvent PRIMME_GLOBAL_SUM;
   PetscLogEventBegin(PRIMME_GLOBAL_SUM,0,0,0,0);
   *ierr = MPI_Allreduce(sendBuf, recvBuf, *count, MPI_DOUBLE, MPI_SUM, communicator);
   PetscLogEventEnd(PRIMME_GLOBAL_SUM,0,0,0,0);
#else
   *ierr = MPI_Allreduce(sendBuf, recvBuf, *count, MPI_DOUBLE, MPI_SUM, communicator);
#endif
}


#endif /* USE_MPI */

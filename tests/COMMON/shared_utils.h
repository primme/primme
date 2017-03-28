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
 * File: shared_utils.h
 * 
 * Purpose - Definitions in shared_utils.c used by the driver.
 * 
 ******************************************************************************/

#ifndef DRIVER_SHARED_H
#define DRIVER_SHARED_H

#include "primme.h"
#include "num.h"

#define TRUE  1
#define FALSE 0
#define ZERO 0.0L

typedef enum {
   driver_default,
   driver_native,
   driver_petsc,
   driver_parasails,
   driver_rsb
} driver_mat;

typedef enum {
   driver_noprecond,    /* no preconditioning */
   driver_jacobi,       /* K=Diag(A-shift),   shift provided once by user */
   driver_jacobi_i,     /* Diag(A-shift_i), shifts provided by primme every step */
   driver_ilut,         /* ILUT(A-shift)  , shift provided once by user */
   driver_normal,       /* precond based on A*A, only for SVD */
   driver_bjacobi       /* block jacobi */
} driver_prec;

typedef struct driver_params {

   char outputFileName[512];
   char partId[256];
   char testId[256];
   char partDir[1024];
   char matrixFileName[1024];
   char initialGuessesFileName[1024];
   char saveXFileName[1024];
   double initialGuessesPert;
   char checkXFileName[1024];
   int checkInterface;

   driver_mat matrixChoice;

   int weightedPart;

   /* Preconditioning paramaters for various preconditioners */
   driver_prec PrecChoice;
   int isymm;
   int level;
   double threshold;
   double filter;
   double shift;
   
} driver_params;


int read_solver_params(char *configFileName, char *outputFileName,
                primme_params *primme, const char* primmeprefix,
                primme_preset_method *method, const char* methodstr);
int read_solver_params_svds(char *configFileName, char *outputFileName,
                primme_svds_params *primme, const char* primmeprefix,
                primme_svds_preset_method *method, const char* methodstr,
                primme_preset_method *primme_method,
                primme_preset_method *primme_method0);
int read_driver_params(char *configFileName, driver_params *driver);
void driver_display_params(driver_params driver, FILE *outputFile);
void driver_display_method(primme_preset_method method, const char *methodstr, FILE *outputFile);
void driver_display_methodsvd(primme_svds_preset_method method, const char *methodstr, FILE *outputFile);
#ifdef USE_MPI
#include <mpi.h>
void par_GlobalSumDouble(void *sendBuf, void *recvBuf, int *count, 
                         primme_params *primme, int *ierr);
void par_GlobalSumDoubleSvds(void *sendBuf, void *recvBuf, int *count, 
                         primme_svds_params *primme, int *ierr);
void broadCast_svds(primme_svds_params *primme_svds, primme_svds_preset_method *method,
   primme_preset_method *primmemethod, primme_preset_method *primmemethod0,
   driver_params *driver, int master, MPI_Comm comm);
void broadCast(primme_params *primme, primme_preset_method *method, 
   driver_params *driver, int master, MPI_Comm comm);

#endif

#endif /* DRIVER_SHARED_H */

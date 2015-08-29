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
 * File: errors.c
 *
 * Purpose - Maintains & displays a stack trace for reporting errors.
 *
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "primme.h"
#include "errors_private.h"


/*******************************************************************************
 * Subroutine primme_PushErrorMessage - Pushes an error message onto the stack
 *   trace.
 * 
 * Parameters
 * ----------
 * callingFunction  Integer code for the calling routine.  This is the routine
 *                  pushing the error message onto the stack.
 *
 * failedFunction   Integer code for the function that returned the error.
 *
 * errorCode        The error code returned by the failing function
 *
 * fileName         The name of the source file the error occured in
 *
 * lineNumber       The line number within the source file where the error
 *                  was trapped
 *
 * Input/Output parameters
 * -----------------------
 * primme_params      Pointer to the primme parameter structure
 *
 ******************************************************************************/

void primme_PushErrorMessage(const primme_function callingFunction, 
     const primme_function failedFunction, const int errorCode, 
     const char *fileName, const int lineNumber, primme_params *primme) {

   stackTraceNode *newNode;

   /* Allocate and initialize a new stack trace node */

   newNode =(stackTraceNode *)primme_calloc(1,sizeof(stackTraceNode),"newNode");
   newNode->callingFunction = callingFunction;
   newNode->failedFunction = failedFunction;
   newNode->errorCode = errorCode;
   newNode->lineNumber = lineNumber;
   strncpy(newNode->fileName, fileName, PRIMME_MAX_NAME_LENGTH);

   /* Push the new node onto the stack */

   if (primme->stackTrace == NULL) {
      primme->stackTrace = newNode;
      newNode->nextNode = NULL;
   }
   else {
      newNode->nextNode = primme->stackTrace;
      primme->stackTrace = newNode;
   }

}


/*******************************************************************************
 * Subroutine primme_PrintStackTrace - This subroutine prints out the stack
 *   trace.
 *
 * Parameters
 * ----------
 * primme_params - primme parameter structure
 *
 ******************************************************************************/
 
void primme_PrintStackTrace(const primme_params primme) {

   char callingFunctionName[PRIMME_MAX_NAME_LENGTH];
   char failingFunctionName[PRIMME_MAX_NAME_LENGTH];

   stackTraceNode *nextNode;
   
   nextNode = primme.stackTrace;
   if (nextNode == NULL) {
      fprintf(stderr, "\n    Successful return\n\n");
   }
   else {

      fprintf(stderr, "\
================================================================================\n                             primme Stack Trace\n\
--------------------------------------------------------------------------------\n");

      /* For each node in the stack, print the error message. */

      while(nextNode != NULL) {
         convertToString(nextNode->callingFunction, callingFunctionName);
         convertToString(nextNode->failedFunction, failingFunctionName);
         fprintf(stderr, 
         "Calling Function: %s File: %s Line: %d Function '%s' returned: %d\n", 
            callingFunctionName, nextNode->fileName, nextNode->lineNumber, 
            failingFunctionName, nextNode->errorCode);
         nextNode = nextNode->nextNode;
      }

   fprintf(stderr, "\
============================== End Stack Trace =================================\n\n");

   }

}
      

/*******************************************************************************
 * Subroutine primme_DeleteStackTrace - This subroutine deletes all the nodes
 *    in the stack.
 *
 * Input/Output parameters
 * -----------------------
 * primme_params - primme parameter structure
 *
 ******************************************************************************/
 
void primme_DeleteStackTrace(primme_params *primme) {

   stackTraceNode *prevNode, *nextNode;

   nextNode = primme->stackTrace;

   while (nextNode != NULL) {
      prevNode = nextNode;
      nextNode = nextNode->nextNode;
      free(prevNode);
   }

   primme->stackTrace = NULL;

}
      
      
/*******************************************************************************
 * Subroutine convertToString - This subroutine accepts an integer code
 *   representing a function.  The code is then converted into the 
 *   corresponding string containing the function name.
 *
 * Parameters
 * ----------
 * functionCode - The code for the function
 *
 * functionName - The name of the function 
 *
 ******************************************************************************/

static void convertToString(primme_function func, char *functionName) {

   switch (func) {
      case Primme_dprimme:
         strcpy(functionName, "dprimme");
         break;
      case Primme_zprimme:
         strcpy(functionName, "zprimme");
         break;
      case Primme_main_iter:
         strcpy(functionName, "main_iter");
         break;
      case Primme_allocate_workspace:
         strcpy(functionName, "allocate_workspace");
         break;
      case Primme_check_input:
         strcpy(functionName, "check_input");
         break;
      case Primme_init_basis:
         strcpy(functionName, "init_basis");
         break;
      case Primme_init_block_krylov:
         strcpy(functionName, "init_block_krylov");
         break;
      case Primme_init_krylov:
         strcpy(functionName, "init_krylov");
         break;
      case Primme_ortho:
         strcpy(functionName, "ortho");
         break;
      case Primme_solve_h:
         strcpy(functionName, "solve_H");
         break;
      case Primme_restart:
         strcpy(functionName, "restart");
         break;
      case Primme_restart_h:
         strcpy(functionName, "restart_H");
         break;
      case Primme_insert_submatrix:
         strcpy(functionName, "insert_submatrix");
         break;
      case Primme_lock_vectors:
         strcpy(functionName, "lock_vectors");
         break;
      case Primme_num_dsyev:
         strcpy(functionName, "Num_dsyev");
         break;
      case Primme_num_zheev:
         strcpy(functionName, "Num_zheev");
         break;
      case Primme_num_dspev:
         strcpy(functionName, "Num_dspev");
         break;
      case Primme_num_zhpev:
         strcpy(functionName, "Num_zhpev");
         break;
      case Primme_ududecompose:
         strcpy(functionName, "UDUDecompose");
         break;
      case Primme_udusolve:
         strcpy(functionName, "UDUSolve");
         break;
      case Primme_apply_projected_preconditioner:
         strcpy(functionName, "apply_projected_preconditioner");
         break;
      case Primme_apply_skew_projector:
         strcpy(functionName, "apply_skew_projector");
         break;
      case Primme_inner_solve:
         strcpy(functionName, "inner_solve");
         break;
      case Primme_solve_correction:
         strcpy(functionName, "solve_correction");
         break;
      case Primme_fopen:
         strcpy(functionName, "fopen");
         break;
      case Primme_malloc:
         strcpy(functionName, "malloc");
         break;
   }

}


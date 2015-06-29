/* ssrcsr.f -- translated by f2c (version 19970219).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
#include "f2c.h"

 ssrcsr is part of the SPARSKIT library by Y. Saad at University of Minnesota
 For the full SPARSKIT library visit:
 	                  www.cs.umn.edu/~saad/
 
*/

int ssrcsr(int *job, int *value2, int *nrow, double *a, int *ja, int *ia, 
   int *nzmax, double *ao, int *jao, int *iao, int *indu, int *iwk, int *ierr)
{
    /* System generated locals */
    int i__1, i__2, i__3;

    /* Local variables */
    static int ipos, i__, j, k, klast, kosav, ko, kfirst;
    static double tmp;
    static int nnz;

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/* -----------------------------------------------------------------------
 */
/*     Symmetric Sparse Row to Compressed Sparse Row format */
/* -----------------------------------------------------------------------
 */
/*     This subroutine converts a given matrix in SSR format to regular */
/*     CSR format by computing Ao = A + A' - diag(A), where A' is A */
/*     transpose. */

/*     Typically this routine is used to expand the SSR matrix of */
/*     Harwell Boeing matrices, or to obtain a symmetrized graph of */
/*     unsymmetric matrices. */

/*     This routine is inplace, i.e., (Ao,jao,iao) may be same as */
/*     (a,ja,ia). */

/*     It is possible to input an arbitrary CSR matrix to this routine, */
/*     since there is no syntactical difference between CSR and SSR */
/*     format. It also removes duplicate entries and perform a partial */
/*     ordering. The output matrix has an order of lower half, main */
/*     diagonal and upper half after the partial ordering. */
/* -----------------------------------------------------------------------
 */
/* on entry: */
/* --------- */

/* job   = options */
/*         0 -- duplicate entries are not removed. If the input matrix is 
*/
/*             SSR (not an arbitary CSR) matrix, no duplicate entry should
*/
/*              arise from this routine. */
/*         1 -- eliminate duplicate entries, zero entries. */
/*         2 -- eliminate duplicate entries and perform partial ordering. 
*/
/*         3 -- eliminate duplicate entries, sort the entries in the */
/*              increasing order of clumn indices. */

/* value2= will the values of A be copied? */
/*         0 -- only expand the graph (a, ao are not touched) */
/*         1 -- expand the matrix with the values. */

/* nrow  = column dimension of inout matrix */
/* a, */
/* ia, */
/* ja    = matrix in compressed sparse row format. */

/* nzmax = size of arrays ao and jao. SSRCSR will abort if the storage */
/*          provided in ao, jao is not sufficient to store A. See ierr. */

/* on return: */
/* ---------- */
/* ao, jao, iao */
/*       = output matrix in compressed sparse row format. The resulting */
/*         matrix is symmetric and is equal to A+A'-D. ao, jao, iao, */
/*         can be the same as a, ja, ia in the calling sequence. */

/* indu  = integer array of length nrow. INDU will contain pointers */
/*         to the beginning of upper traigular part if job > 1. */
/*         Otherwise it is also used as a work array (size nrow). */

/* iwk   = integer work space (size nrow+1). */

/* ierr  = integer. Serving as error message. If the length of the arrays 
*/
/*         ao, jao exceeds nzmax, ierr returns the minimum value */
/*         needed for nzmax. otherwise ierr=0 (normal return). */

/* -----------------------------------------------------------------------
 */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Executable Statements .. */
    /* Parameter adjustments */
    --iwk;
    --indu;
    --iao;
    --ia;
    --a;
    --ja;
    --jao;
    --ao;

    /* Function Body */
    *ierr = 0;
    i__1 = *nrow;
    for (i__ = 1; i__ <= i__1; ++i__) {
	indu[i__] = 0;
	iwk[i__] = 0;
/* L10: */
    }
    iwk[*nrow + 1] = 0;

/*     .. compute number of elements in each row of (A'-D) */
/*     put result in iwk(i+1)  for row i. */

    i__1 = *nrow;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = ia[i__ + 1] - 1;
	for (k = ia[i__]; k <= i__2; ++k) {
	    j = ja[k];
	    if (j != i__) {
		++iwk[j + 1];
	    }
/* L20: */
	}
/* L30: */
    }

/*     .. find addresses of first elements of ouput matrix. result in iwk 
*/

    iwk[1] = 1;
    i__1 = *nrow;
    for (i__ = 1; i__ <= i__1; ++i__) {
	indu[i__] = iwk[i__] + ia[i__ + 1] - ia[i__];
	iwk[i__ + 1] += indu[i__];
	--indu[i__];
/* L40: */
    }
/* .....Have we been given enough storage in ao, jao ? */
    nnz = iwk[*nrow + 1] - 1;
    if (nnz > *nzmax) {
	*ierr = nnz;
	return 0;
    }

/*     .. copy the existing matrix (backwards). */

    kosav = iwk[*nrow + 1];
    for (i__ = *nrow; i__ >= 1; --i__) {
	klast = ia[i__ + 1] - 1;
	kfirst = ia[i__];
	iao[i__ + 1] = kosav;
	kosav = iwk[i__];
	ko = iwk[i__] - kfirst;
	iwk[i__] = ko + klast + 1;
	i__1 = kfirst;
	for (k = klast; k >= i__1; --k) {
	    if (*value2 != 0) {
		ao[k + ko] = a[k];
	    }
	    jao[k + ko] = ja[k];
/* L50: */
	}
/* L60: */
    }
    iao[1] = 1;

/*     now copy (A'-D). Go through the structure of ao, jao, iao */
/*     that has already been copied. iwk(i) is the address */
/*     of the next free location in row i for ao, jao. */

    i__1 = *nrow;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = indu[i__];
	for (k = iao[i__]; k <= i__2; ++k) {
	    j = jao[k];
	    if (j != i__) {
		ipos = iwk[j];
		if (*value2 != 0) {
		    ao[ipos] = ao[k];
		}
		jao[ipos] = i__;
		iwk[j] = ipos + 1;
	    }
/* L70: */
	}
/* L80: */
    }
    if (*job <= 0) {
	return 0;
    }

/*     .. eliminate duplicate entries -- */
/*     array INDU is used as marker for existing indices, it is also the 
*/
/*     location of the entry. */
/*     IWK is used to stored the old IAO array. */
/*     matrix is copied to squeeze out the space taken by the duplicated 
*/
/*     entries. */

    i__1 = *nrow;
    for (i__ = 1; i__ <= i__1; ++i__) {
	indu[i__] = 0;
	iwk[i__] = iao[i__];
/* L90: */
    }
    iwk[*nrow + 1] = iao[*nrow + 1];
    k = 1;
    i__1 = *nrow;
    for (i__ = 1; i__ <= i__1; ++i__) {
	iao[i__] = k;
	ipos = iwk[i__];
	klast = iwk[i__ + 1];
L100:
	if (ipos < klast) {
	    j = jao[ipos];
	    if (indu[j] == 0) {
/*     .. new entry .. */
		if (*value2 != 0) {
		    if (ao[ipos] != 0.) {
			indu[j] = k;
			jao[k] = jao[ipos];
			ao[k] = ao[ipos];
			++k;
		    }
		} else {
		    indu[j] = k;
		    jao[k] = jao[ipos];
		    ++k;
		}
	    } else if (*value2 != 0) {
/*     .. duplicate entry .. */
		ao[indu[j]] += ao[ipos];
	    }
	    ++ipos;
	    goto L100;
	}
/*     .. remove marks before working on the next row .. */
	i__2 = k - 1;
	for (ipos = iao[i__]; ipos <= i__2; ++ipos) {
	    indu[jao[ipos]] = 0;
/* L110: */
	}
/* L120: */
    }
    iao[*nrow + 1] = k;
    if (*job <= 1) {
	return 0;
    }

/*     .. partial ordering .. */
/*     split the matrix into strict upper/lower triangular */
/*     parts, INDU points to the the beginning of the strict upper part. 
*/

    i__1 = *nrow;
    for (i__ = 1; i__ <= i__1; ++i__) {
	klast = iao[i__ + 1] - 1;
	kfirst = iao[i__];
L130:
	if (klast > kfirst) {
	    if (jao[klast] < i__ && jao[kfirst] >= i__) {
/*     .. swap klast with kfirst .. */
		j = jao[klast];
		jao[klast] = jao[kfirst];
		jao[kfirst] = j;
		if (*value2 != 0) {
		    tmp = ao[klast];
		    ao[klast] = ao[kfirst];
		    ao[kfirst] = tmp;
		}
	    }
	    if (jao[klast] >= i__) {
		--klast;
	    }
	    if (jao[kfirst] < i__) {
		++kfirst;
	    }
	    goto L130;
	}

	if (jao[klast] < i__) {
	    indu[i__] = klast + 1;
	} else {
	    indu[i__] = klast;
	}
/* L140: */
    }
    if (*job <= 2) {
	return 0;
    }

/*     .. order the entries according to column indices */
/*     bubble-sort is used */

    i__1 = *nrow;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = indu[i__] - 1;
	for (ipos = iao[i__]; ipos <= i__2; ++ipos) {
	    i__3 = ipos + 1;
	    for (j = indu[i__] - 1; j >= i__3; --j) {
		k = j - 1;
		if (jao[k] > jao[j]) {
		    ko = jao[k];
		    jao[k] = jao[j];
		    jao[j] = ko;
		    if (*value2 != 0) {
			tmp = ao[k];
			ao[k] = ao[j];
			ao[j] = tmp;
		    }
		}
/* L150: */
	    }
/* L160: */
	}
	i__2 = iao[i__ + 1] - 1;
	for (ipos = indu[i__]; ipos <= i__2; ++ipos) {
	    i__3 = ipos + 1;
	    for (j = iao[i__ + 1] - 1; j >= i__3; --j) {
		k = j - 1;
		if (jao[k] > jao[j]) {
		    ko = jao[k];
		    jao[k] = jao[j];
		    jao[j] = ko;
		    if (*value2 != 0) {
			tmp = ao[k];
			ao[k] = ao[j];
			ao[j] = tmp;
		    }
		}
/* L170: */
	    }
/* L180: */
	}
/* L190: */
    }

    return 0;
/* ---- end of ssrcsr ----------------------------------------------------
 */
/* -----------------------------------------------------------------------
 */
} /* ssrcsr_ */




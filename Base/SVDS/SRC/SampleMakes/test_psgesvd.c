#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include "mpi.h"
#include <complex.h>
#include "blas.h"
#include "blacs.h"
#include "lapack.h"
#include "scalapack.h"

static int max( int a, int b ){
        if (a>b) return(a); else return(b);
}
static int min( int a, int b ){
        if (a<b) return(a); else return(b);
}

extern float verif_orthogonality(int m, int n, float *U, int iu, int ju, int *descU);
extern float verif_representativity(int m, int n, float *A, int ia, int ja, int *descA,
                                                     float *U, int iu, int ju, int *descU,
                                                     float *VT, int ivt, int jvt, int *descVT,
                                                     float *S);
extern float verif_repres_NV(int m, int n, float *A, int ia, int ja, int *descA,
                                              float *VT, int ivt, int jvt, int *descVT,
                                              float *S);
extern float verif_repres_VN(int m, int n, float *A, int ia, int ja, int *descA,
                                              float *U, int iu, int ju, int *descU,
                                              float *S);
extern int driver_psgesvd(char jobU, char jobVT, int m, int n, float *A, int ia, int ja, int *descA,
			float *S_NN, float *U_NN, int iu, int ju, int *descU, float *VT_NN, int ivt, int jvt, int *descVT,
		       	double *MPIelapsedNN);


int main(int argc, char **argv) {
	int iam, nprocs;
	int myrank_mpi, nprocs_mpi;
	int ictxt, nprow, npcol, myrow, mycol;
	int nb, m, n;
	int mpA, nqA, mpU, nqU, mpVT, nqVT;
	int i, j, k, itemp, min_mn;
	int descA[9], descU[9], descVT[9];
	float *A=NULL;
	int info, infoNN, infoVV, infoNV, infoVN;
	float *U_NN=NULL,  *U_VV=NULL,  *U_NV=NULL,  *U_VN=NULL;
	float *VT_NN=NULL, *VT_VV=NULL, *VT_NV=NULL, *VT_VN=NULL;
	float *S_NN=NULL,  *S_VV=NULL, *S_NV=NULL, *S_VN=NULL;
	float *S_res_NN=NULL;
        float orthU_VV, residF, orthVT_VV;
        float orthU_VN, orthVT_NV;
	float  residS_NN, eps;
	float  res_repres_NV, res_repres_VN;
/**/
	int izero=0,ione=1;
	float rtmone=-1.0e+00;
/**/
	double MPIelapsedVV, MPIelapsedNN, MPIelapsedVN, MPIelapsedNV;
	char jobU, jobVT;
	int nbfailure=0, nbtestcase=0,inputfromfile, nbhetereogeneity=0;
	float threshold=100e+00;
	char buf[1024];
	FILE *fd;	
	char *c;
	char *t_jobU, *t_jobVT;
	int *t_m, *t_n, *t_nb, *t_nprow, *t_npcol;
	int nb_expe, expe;
	char hetereogeneityVV, hetereogeneityNN, hetereogeneityVN, hetereogeneityNV;
	int iseed[4], idist;
/**/
	MPI_Init( &argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank_mpi);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs_mpi);
/**/
	m = 100; n = 100; nprow = 1; npcol = 1; nb = 64; jobU='A'; jobVT='A'; inputfromfile = 0;
	for( i = 1; i < argc; i++ ) {
		if( strcmp( argv[i], "-f" ) == 0 ) {
			inputfromfile = 1;
		}
		if( strcmp( argv[i], "-jobvt" ) == 0 ) {
			if (i+1<argc) {
				if( strcmp( argv[i+1], "V" ) == 0 ){ jobVT = 'V'; i++; }
				else if( strcmp( argv[i+1], "N" ) == 0 ){ jobVT = 'N'; i++; }
				else if( strcmp( argv[i+1], "A" ) == 0 ){ jobVT = 'A'; i++; }
				else printf(" ** warning: jobvt should be set to V, N or A in the command line ** \n");
			}
			else	
				printf(" ** warning: jobvt should be set to V, N or A in the command line ** \n");
		}
		if( strcmp( argv[i], "-jobu" ) == 0 ) {
			if (i+1<argc) {
				if( strcmp( argv[i+1], "V" ) == 0 ){ jobU = 'V'; i++; }
				else if( strcmp( argv[i+1], "N" ) == 0 ){ jobU = 'N'; i++; }
				else if( strcmp( argv[i+1], "A" ) == 0 ){ jobU = 'A'; i++; }
				else printf(" ** warning: jobu should be set to V, N or A in the command line ** \n");
			}
			else	
				printf(" ** warning: jobu should be set to V, N or A in the command line ** \n");
		}
		if( strcmp( argv[i], "-m" ) == 0 ) {
			m      = atoi(argv[i+1]);
			i++;
		}
		if( strcmp( argv[i], "-n" ) == 0 ) {
			n      = atoi(argv[i+1]);
			i++;
		}
		if( strcmp( argv[i], "-p" ) == 0 ) {
			nprow  = atoi(argv[i+1]);
			i++;
		}
		if( strcmp( argv[i], "-q" ) == 0 ) {
			npcol  = atoi(argv[i+1]);
			i++;
		}
		if( strcmp( argv[i], "-nb" ) == 0 ) {
			nb     = atoi(argv[i+1]);
			i++;
		}
	}
/**/
	if (inputfromfile){
		nb_expe = 0;
		fd = fopen("svd.dat", "r");
		if (fd == NULL) { printf("File failed to open svd.dat from processor mpirank(%d/%d): \n",myrank_mpi,nprocs_mpi); exit(-1); }
		do {	
			c = fgets(buf, 1024, fd);  /* get one line from the file */
			if (c != NULL)
				if (c[0] != '#')
					nb_expe++;
		} while (c != NULL);              /* repeat until NULL          */
		fclose(fd);
		t_jobU  = (char *)calloc(nb_expe,sizeof(char)) ;
		t_jobVT = (char *)calloc(nb_expe,sizeof(char)) ;
		t_m     = (int  *)calloc(nb_expe,sizeof(int )) ;
		t_n     = (int  *)calloc(nb_expe,sizeof(int )) ;
		t_nb    = (int  *)calloc(nb_expe,sizeof(int )) ;
		t_nprow = (int  *)calloc(nb_expe,sizeof(int )) ;
		t_npcol = (int  *)calloc(nb_expe,sizeof(int )) ;
		fd = fopen("svd.dat", "r");
		expe=0;
		do {	
			c = fgets(buf, 1024, fd);  /* get one line from the file */
			if (c != NULL)
				if (c[0] != '#'){
					//printf("NBEXPE = %d\n",expe);
					sscanf(c,"%c %c %d %d %d %d %d",
						&(t_jobU[expe]),&(t_jobVT[expe]),&(t_m[expe]),&(t_n[expe]),
						&(t_nb[expe]),(&t_nprow[expe]),&(t_npcol[expe]));
					expe++;
				}
		} while (c != NULL);              /* repeat until NULL          */
		fclose(fd);
	}
	else {
		nb_expe = 1;
		t_jobU  = (char *)calloc(nb_expe,sizeof(char)) ;
		t_jobVT = (char *)calloc(nb_expe,sizeof(char)) ;
		t_m     = (int  *)calloc(nb_expe,sizeof(int )) ;
		t_n     = (int  *)calloc(nb_expe,sizeof(int )) ;
		t_nb    = (int  *)calloc(nb_expe,sizeof(int )) ;
		t_nprow = (int  *)calloc(nb_expe,sizeof(int )) ;
		t_npcol = (int  *)calloc(nb_expe,sizeof(int )) ;
		t_jobU[0]  = jobU;
		t_jobVT[0] = jobVT;
		t_m[0]     = m;
		t_n[0]     = n;
		t_nb[0]    = nb;
		t_nprow[0] = nprow;
		t_npcol[0] = npcol;
	}

        if (myrank_mpi==0){
		printf("\n");
		printf("--------------------------------------------------------------------------------------------------------------------\n");
				printf("                            Testing psgsevd -- float precision SVD ScaLAPACK routine                \n");
		printf("jobU jobVT    m     n     nb   p   q   || info   heter   resid     orthU    orthVT   |SNN-SVV|    time(s)   cond(A) \n");
		printf("--------------------------------------------------------------------------------------------------------------------\n");
	}
/**/
	for (expe = 0; expe<nb_expe; expe++){

	jobU  = t_jobU[expe]  ; 
	jobVT = t_jobVT[expe] ; 
	m     = t_m[expe]     ; 
	n     = t_n[expe]     ; 
	nb    = t_nb[expe]    ; 
	nprow = t_nprow[expe] ; 
	npcol = t_npcol[expe] ; 

	if (nb>n)
		nb = n;
	if (nprow*npcol>nprocs_mpi){
		if (myrank_mpi==0)
			printf(" **** ERROR : we do not have enough processes available to make a p-by-q process grid ***\n");
			printf(" **** Bye-bye                                                                         ***\n");
		MPI_Finalize(); exit(1);
	}
/**/
	Cblacs_pinfo( &iam, &nprocs ) ;
	Cblacs_get( -1, 0, &ictxt );
	Cblacs_gridinit( &ictxt, "Row", nprow, npcol );
	Cblacs_gridinfo( ictxt, &nprow, &npcol, &myrow, &mycol );
/**/
	min_mn = min(m,n);
/**/
 	//if (iam==0)
 		//printf("\tm=%d\tn = %d\t\t(%d,%d)\t%dx%d\n",m,n,nprow,npcol,nb,nb);
 	//printf("Hello World, I am proc %d over %d for MPI, proc %d over %d for BLACS in position (%d,%d) in the process grid\n", 
 	 		//myrank_mpi,nprocs_mpi,iam,nprocs,myrow,mycol);
/*
*
*     Work only the process in the process grid
*
*/
	//if ((myrow < nprow)&(mycol < npcol)){
	if ((myrow>-1)&(mycol>-1)&(myrow<nprow)&(mycol<npcol)){

/*
*
*     Compute the size of the local matrices (thanks to numroc)
*
*/ 
		mpA    = numroc_( &m     , &nb, &myrow, &izero, &nprow );
		nqA    = numroc_( &n     , &nb, &mycol, &izero, &npcol );
		mpU    = numroc_( &m     , &nb, &myrow, &izero, &nprow );
		nqU    = numroc_( &min_mn, &nb, &mycol, &izero, &npcol );
		mpVT   = numroc_( &min_mn, &nb, &myrow, &izero, &nprow );
		nqVT   = numroc_( &n     , &nb, &mycol, &izero, &npcol );
/*
*
*     Allocate and fill the matrices A and B
*
*/ 
		A = (float *)calloc(mpA*nqA,sizeof(float)) ;
		if (A==NULL){ printf("error of memory allocation A on proc %dx%d\n",myrow,mycol); exit(0); }
/**/		
//		seed = iam*(mpA*nqA*2); srand(seed);
		idist = 2;
		iseed[0] = mpA%4096;
		iseed[1] = iam%4096;
		iseed[2] = nqA%4096;
		iseed[3] = 23;
/**/		
		k = 0;
		for (i = 0; i < mpA; i++) {
			for (j = 0; j < nqA; j++) {
				slarnv_( &idist, iseed, &ione, &(A[k]) );
				k++;	
			}
		}
/*
*
*     Initialize the array descriptor for the distributed matrices xA, U and VT
*
*/ 
		itemp = max( 1, mpA );
		descinit_( descA,  &m, &n, &nb, &nb, &izero, &izero, &ictxt, &itemp, &info );
		itemp = max( 1, mpA );
		descinit_( descU,  &m, &min_mn, &nb, &nb, &izero, &izero, &ictxt, &itemp, &info );
		itemp = max( 1, mpVT );
		descinit_( descVT, &min_mn, &n, &nb, &nb, &izero, &izero, &ictxt, &itemp, &info );
/**/
      		eps = pslamch_( &ictxt, "Epsilon" );
/**/
		if ( ((jobU=='V')&(jobVT=='N')) ||(jobU == 'A' )||(jobVT=='A')){
		nbtestcase++;	
		U_VN = (float *)calloc(mpU*nqU,sizeof(float)) ;
		if (U_VN==NULL){ printf("error of memory allocation U_VN on proc %dx%d\n",myrow,mycol); exit(0); }
		S_VN = (float *)calloc(min_mn,sizeof(float)) ;
		if (S_VN==NULL){ printf("error of memory allocation S_VN on proc %dx%d\n",myrow,mycol); exit(0); }
		infoVN = driver_psgesvd( 'V', 'N', m, n, A, 1, 1, descA,
			S_VN, U_VN, 1, 1, descU, VT_VN, 1, 1, descVT,
		       	&MPIelapsedVN);
		orthU_VN  = verif_orthogonality(m,min_mn,U_VN , 1, 1, descU);
		res_repres_VN = verif_repres_VN( m, n, A, 1, 1, descA, U_VN, 1, 1, descU, S_VN);
		if (infoVN==min_mn+1) hetereogeneityVN = 'H'; else hetereogeneityVN = 'N';
		if ( iam==0 )
			printf(" V    N   %6d %6d  %3d  %3d %3d  ||  %3d     %c   %7.1e   %7.1e                        %8.2f    %7.1e\n",
				m,n,nb,nprow,npcol,infoVN,hetereogeneityVN,res_repres_VN/(S_VN[0]/S_VN[min_mn-1]),
				orthU_VN,MPIelapsedVN,S_VN[0]/S_VN[min_mn-1]);
		if (infoVN==min_mn+1) nbhetereogeneity++ ;
		else if ((res_repres_VN/eps/(S_VN[0]/S_VN[min_mn-1])>threshold)||(orthU_VN/eps>threshold)||(infoVN!=0)) nbfailure++;
		}
/**/
		if (((jobU=='N')&(jobVT=='V'))||(jobU == 'A' )||(jobVT=='A')){
		nbtestcase++;	
		VT_NV = (float *)calloc(mpVT*nqVT,sizeof(float)) ;
		if (VT_NV==NULL){ printf("error of memory allocation VT_NV on proc %dx%d\n",myrow,mycol); exit(0); }
		S_NV = (float *)calloc(min_mn,sizeof(float)) ;
		if (S_NV==NULL){ printf("error of memory allocation S_NV on proc %dx%d\n",myrow,mycol); exit(0); }
		infoNV = driver_psgesvd( 'N', 'V', m, n, A, 1, 1, descA,
			S_NV, U_NV, 1, 1, descU, VT_NV, 1, 1, descVT,
		       	&MPIelapsedNV);
		orthVT_NV = verif_orthogonality(min_mn,n,VT_NV, 1, 1, descVT);
		res_repres_NV = verif_repres_NV( m, n, A, 1, 1, descA, VT_NV, 1, 1, descVT, S_NV);
		if (infoNV==min_mn+1) hetereogeneityNV = 'H'; else hetereogeneityNV = 'N';
		if ( iam==0 )
			printf(" N    V   %6d %6d  %3d  %3d %3d  ||  %3d     %c   %7.1e             %7.1e              %8.2f    %7.1e\n",
				m,n,nb,nprow,npcol,infoNV,hetereogeneityNV,res_repres_NV/(S_NV[0]/S_NV[min_mn-1]),
				orthVT_NV,MPIelapsedNV,S_NV[0]/S_NV[min_mn-1]);
		if (infoNV==min_mn+1) nbhetereogeneity++ ;
		else if ((res_repres_NV/eps/(S_NV[0]/S_NV[min_mn-1])>threshold)||(orthVT_NV/eps>threshold)||(infoNV!=0)) nbfailure++;
		}
/**/
		if ( ((jobU=='N')&(jobVT=='N')) || ((jobU=='V')&(jobVT=='V')) || (jobU == 'A' ) || (jobVT=='A') ) {
		nbtestcase++;	
		U_VV = (float *)calloc(mpU*nqU,sizeof(float)) ;
		if (U_VV==NULL){ printf("error of memory allocation U_VV on proc %dx%d\n",myrow,mycol); exit(0); }
		VT_VV = (float *)calloc(mpVT*nqVT,sizeof(float)) ;
		if (VT_VV==NULL){ printf("error of memory allocation VT_VV on proc %dx%d\n",myrow,mycol); exit(0); }
		S_VV = (float *)calloc(min_mn,sizeof(float)) ;
		if (S_VV==NULL){ printf("error of memory allocation S_VV on proc %dx%d\n",myrow,mycol); exit(0); }
		infoVV = driver_psgesvd( 'V', 'V', m, n, A, 1, 1, descA,
			S_VV, U_VV, 1, 1, descU, VT_VV, 1, 1, descVT,
		       	&MPIelapsedVV);
		orthU_VV  = verif_orthogonality(m,min_mn,U_VV , 1, 1, descU);
		orthVT_VV = verif_orthogonality(min_mn,n,VT_VV, 1, 1, descVT);
		residF =  verif_representativity( m, n,     A, 1, 1, descA,
                                                         U_VV, 1, 1, descU,
                                                        VT_VV, 1, 1, descVT,
                                                         S_VV);
		if (infoVV==min_mn+1) hetereogeneityVV = 'H'; else hetereogeneityVV = 'N';
		if ( iam==0 )
			printf(" V    V   %6d %6d  %3d  %3d %3d  ||  %3d     %c   %7.1e   %7.1e   %7.1e              %8.2f    %7.1e\n",
				m,n,nb,nprow,npcol,infoVV,hetereogeneityVV,residF,orthU_VV,orthVT_VV,MPIelapsedVV,S_VV[0]/S_VV[min_mn-1]);
		if (infoVV==min_mn+1) nbhetereogeneity++ ;
		else if ((residF/eps>threshold)||(orthU_VV/eps>threshold)||(orthVT_VV/eps>threshold)||(infoVV!=0)) nbfailure++;
		}
/**/
		if (((jobU=='N')&(jobVT=='N'))||(jobU == 'A' )||(jobVT=='A')){
		nbtestcase++;	
		S_NN = (float *)calloc(min_mn,sizeof(float)) ;
		if (S_NN==NULL){ printf("error of memory allocation S_NN on proc %dx%d\n",myrow,mycol); exit(0); }
		infoNN = driver_psgesvd( 'N', 'N', m, n, A, 1, 1, descA,
			S_NN, U_NN, 1, 1, descU, VT_NN, 1, 1, descVT,
		       	&MPIelapsedNN);
		S_res_NN = (float *)calloc(min_mn,sizeof(float)) ;
		if (S_res_NN==NULL){ printf("error of memory allocation S on proc %dx%d\n",myrow,mycol); exit(0); }
		scopy_(&min_mn,S_VV,&ione,S_res_NN,&ione);
		saxpy_ (&min_mn,&rtmone,S_NN,&ione,S_res_NN,&ione);
		residS_NN = snrm2_(&min_mn,S_res_NN,&ione) / snrm2_(&min_mn,S_VV,&ione);
		free(S_res_NN);
		if (infoNN==min_mn+1) hetereogeneityNN = 'H'; else hetereogeneityNN = 'N';
		if ( iam==0 )
			printf(" N    N   %6d %6d  %3d  %3d %3d  ||  %3d     %c                                  %7.1e   %8.2f    %7.1e\n",
				m,n,nb,nprow,npcol,infoNN,hetereogeneityNN,residS_NN,MPIelapsedNN,S_NN[0]/S_NN[min_mn-1]);
		if (infoNN==min_mn+1) nbhetereogeneity++ ;
		else if ((residS_NN/eps>threshold)||(infoNN!=0)) nbfailure++;
		}
/**/
		if (((jobU=='V')&(jobVT=='N'))||(jobU == 'A' )||(jobVT=='A')){ free(S_VN); free(U_VN); }
		if (((jobU=='N')&(jobVT=='V'))||(jobU == 'A' )||(jobVT=='A')){ free(VT_NV); free(S_NV); }
		if (((jobU=='N')&(jobVT=='N'))||(jobU == 'A' )||(jobVT=='A')){ free(S_NN); }
		if (((jobU=='N')&(jobVT=='N'))||((jobU=='V')&(jobVT=='V'))||(jobU == 'A' )||(jobVT=='A')){ free(U_VV); free(S_VV); free(VT_VV);}
		free(A);
		Cblacs_gridexit( 0 );
	}
/*
*     Print ending messages
*/
	}
	if ( iam==0 ){
		printf("--------------------------------------------------------------------------------------------------------------------\n");
		printf("               [ nbhetereogeneity = %d / %d ]\n",nbhetereogeneity, nbtestcase);
		printf("               [ nbfailure        = %d / %d ]\n",nbfailure, nbtestcase-nbhetereogeneity);
		printf("--------------------------------------------------------------------------------------------------------------------\n");
		printf("\n");
	}
/**/
	free(t_jobU  );
	free(t_jobVT );
	free(t_m     );
	free(t_n     );
	free(t_nb    );
	free(t_nprow );
	free(t_npcol );
	MPI_Finalize();
	exit(0);
}
/**/		
float verif_orthogonality(int m, int n, float *U, int iu, int ju, int *descU){

	float *R=NULL;
	int nprow, npcol, myrow, mycol;
	int mpR, nqR, nb, itemp, descR[9], ictxt, info, min_mn, max_mn;
	int ctxt_ = 1, nb_ = 5;
	int izero = 0, ione = 1;
	float *wwork=NULL;
	float tmone= -1.0e+00,  tpone= +1.0e+00,  tzero= +0.0e+00;
	float orthU;

	min_mn = min(m,n);
	max_mn = max(m,n);
	ictxt = descU[ctxt_];
	nb = descU[nb_];
	Cblacs_gridinfo( ictxt, &nprow, &npcol, &myrow, &mycol );

	mpR    = numroc_( &min_mn, &nb, &myrow, &izero, &nprow );
	nqR    = numroc_( &min_mn, &nb, &mycol, &izero, &npcol );
	R = (float *)calloc(mpR*nqR,sizeof(float)) ;
	if (R==NULL){ printf("error of memory allocation R on proc %dx%d\n",myrow,mycol); exit(0); }
	itemp = max( 1, mpR );
	descinit_( descR,  &min_mn, &min_mn, &nb, &nb, &izero, &izero, &ictxt, &itemp, &info );

	pslaset_( "F", &min_mn, &min_mn, &tzero, &tpone, R, &ione, &ione, descR );
	if (m>n)
		psgemm_( "T", "N", &min_mn, &min_mn, &m, &tpone, U, &iu, &ju, descU, U, 
			&iu, &ju, descU, &tmone, R, &ione, &ione, descR );
	else
		psgemm_( "N", "T", &min_mn, &min_mn, &n, &tpone, U, &iu, &ju, descU, U, 
			&iu, &ju, descU, &tmone, R, &ione, &ione, descR );
	orthU = pslange_( "F", &min_mn, &min_mn, R, &ione, &ione, descR, wwork );
	orthU = orthU / ((float) max_mn);
	free(R);

	return orthU;

}
/**/		
float verif_representativity(int m, int n, float *A, int ia, int ja, int *descA,
                                              float *U, int iu, int ju, int *descU,
                                              float *VT, int ivt, int jvt, int *descVT,
                                              float *S){

	float *Acpy=NULL, *Ucpy=NULL;
	int nprow, npcol, myrow, mycol;
	int min_mn, max_mn, mpA, pcol, localcol, i, nqA;
	int ictxt, nbA, rsrcA, csrcA, nbU, rsrcU, csrcU, mpU, nqU;
	int ctxt_ = 1, nb_ = 5, rsrc_ = 6, csrc_ = 7;
	int izero = 0, ione = 1;
	float *wwork=NULL;
	float tmone= -1.0e+00, tpone= +1.0e+00;
	float residF, AnormF;

	min_mn = min(m,n);
	max_mn = max(m,n);
	ictxt = descA[ctxt_];
	Cblacs_gridinfo( ictxt, &nprow, &npcol, &myrow, &mycol );

	nbA = descA[nb_]; rsrcA = descA[rsrc_] ; csrcA = descA[csrc_] ;
	mpA    = numroc_( &m     , &nbA, &myrow, &rsrcA, &nprow );
	nqA    = numroc_( &n     , &nbA, &mycol, &csrcA, &npcol );
	Acpy = (float *)calloc(mpA*nqA,sizeof(float)) ;
	if (Acpy==NULL){ printf("error of memory allocation Acpy on proc %dx%d\n",myrow,mycol); exit(0); }
      	pslacpy_( "All", &m, &n, A, &ia, &ja, descA, Acpy, &ia, &ja, descA );

	nbU = descU[nb_]; rsrcU = descU[rsrc_] ; csrcU = descU[csrc_] ;
	mpU    = numroc_( &m     , &nbU, &myrow, &rsrcU, &nprow );
	nqU    = numroc_( &min_mn, &nbU, &mycol, &csrcU, &npcol );
	Ucpy = (float *)calloc(mpU*nqU,sizeof(float)) ;
	if (Ucpy==NULL){ printf("error of memory allocation Ucpy on proc %dx%d\n",myrow,mycol); exit(0); }
      	pslacpy_( "All", &m, &min_mn, U, &iu, &ju, descU, Ucpy, &iu, &ju, descU );

	AnormF = pslange_( "F", &m, &n, A, &ia, &ja, descA, wwork);

	for (i=1;i<min_mn+1;i++){
		pcol = indxg2p_( &i, &(descU[5]), &izero, &izero, &npcol );
		localcol = indxg2l_( &i, &(descU[5]), &izero, &izero, &npcol );
		if( mycol==pcol )
			sscal_( &mpA, &(S[i-1]), &(Ucpy[ ( localcol-1 )*descU[8] ]), &ione );
	}
	psgemm_( "N", "N", &m, &n, &min_mn, &tpone, Ucpy, &iu, &ju, descU, VT, &ivt, &jvt, descVT,
	       		&tmone, Acpy, &ia, &ja, descA ); 
	residF = pslange_( "F", &m, &n, Acpy, &ione, &ione, descA, wwork);
	residF = residF/AnormF/((float) max_mn);

	free(Ucpy);
	free(Acpy);

	return residF;
}
/**/		
float verif_repres_NV(int m, int n, float *A, int ia, int ja, int *descA,
                                float *VT, int ivt, int jvt, int *descVT,
                                float *S){

	float *Ucpy=NULL;
	int nprow, npcol, myrow, mycol;
	int min_mn, max_mn, mpA, pcol, localcol, i, nqA;
	int ictxt, nbA, rsrcA, csrcA, mpU, nqU, descUcpy[9], itemp, iucpy, jucpy;
	int ctxt_ = 1, nb_ = 5, rsrc_ = 6, csrc_ = 7;
	int izero = 0, ione = 1, info;
	float tpone= +1.0e+00,  tzero= +0.0e+00;
	float verif_repres_NV, invStemp;

	min_mn = min(m,n);
	max_mn = max(m,n);
	ictxt = descA[ctxt_];
	Cblacs_gridinfo( ictxt, &nprow, &npcol, &myrow, &mycol );

	nbA = descA[nb_]; rsrcA = descA[rsrc_] ; csrcA = descA[csrc_] ;

	mpA    = numroc_( &m     , &nbA, &myrow, &rsrcA, &nprow );
	nqA    = numroc_( &n     , &nbA, &mycol, &csrcA, &npcol );

	itemp = max( 1, mpA );
	descinit_( descUcpy,  &m, &min_mn, &nbA, &nbA, &rsrcA, &csrcA, &ictxt, &itemp, &info );

	iucpy = 1; jucpy = 1;	
	mpU    = numroc_( &m     , &nbA, &myrow, &rsrcA, &nprow );
	nqU    = numroc_( &min_mn, &nbA, &mycol, &csrcA, &npcol );
	Ucpy = (float *)calloc(mpU*nqU,sizeof(float)) ;
	if (Ucpy==NULL){ printf("error of memory allocation Ucpy on proc %dx%d\n",myrow,mycol); exit(0); }

	psgemm_( "N", "T", &m, &min_mn, &n, &tpone, A, &ia, &ja, descA, VT, &ivt, &jvt, descVT,
	       		&tzero, Ucpy, &iucpy, &jucpy, descUcpy ); 

	for (i=1;i<min_mn+1;i++){
		pcol = indxg2p_( &i, &(descUcpy[5]), &izero, &izero, &npcol );
		localcol = indxg2l_( &i, &(descUcpy[5]), &izero, &izero, &npcol );
		invStemp = 1/S[i-1];
		if( mycol==pcol )
			sscal_( &mpA, &invStemp, &(Ucpy[ ( localcol-1 )*descUcpy[8] ]), &ione );
	}

	verif_repres_NV = verif_orthogonality(m,min_mn,Ucpy, iucpy, jucpy, descUcpy);

	free(Ucpy);

	return verif_repres_NV;
}

/**/		
float verif_repres_VN(int m, int n, float *A, int ia, int ja, int *descA,
                                       float *U, int iu, int ju, int *descU,
                                       float *S){

	float *VTcpy=NULL;
	int nprow, npcol, myrow, mycol;
	int min_mn, max_mn, mpA, prow, localcol, i, nqA;
	int ictxt, nbA, rsrcA, csrcA, mpVT, nqVT, descVTcpy[9], itemp, ivtcpy, jvtcpy;
	int ctxt_ = 1, nb_ = 5, rsrc_ = 6, csrc_ = 7;
	int izero = 0, info;
	float tpone= +1.0e+00,  tzero= +0.0e+00;
	float verif_repres_VN, invStemp;

	min_mn = min(m,n);
	max_mn = max(m,n);
	ictxt = descA[ctxt_];
	Cblacs_gridinfo( ictxt, &nprow, &npcol, &myrow, &mycol );

	nbA = descA[nb_]; rsrcA = descA[rsrc_] ; csrcA = descA[csrc_] ;

	mpA    = numroc_( &m     , &nbA, &myrow, &rsrcA, &nprow );
	nqA    = numroc_( &n     , &nbA, &mycol, &csrcA, &npcol );

	mpVT   = numroc_( &min_mn, &nbA, &myrow, &rsrcA, &nprow );
	nqVT   = numroc_( &n     , &nbA, &mycol, &csrcA, &npcol );

	itemp = max( 1, mpVT );
	descinit_( descVTcpy, &min_mn, &n, &nbA, &nbA, &rsrcA, &csrcA, &ictxt, &itemp, &info );

	ivtcpy = 1; jvtcpy = 1;	
	VTcpy = (float *)calloc(mpVT*nqVT,sizeof(float)) ;
	if (VTcpy==NULL){ printf("error of memory allocation VTcpy on proc %dx%d\n",myrow,mycol); exit(0); }

	psgemm_( "T", "N", &min_mn, &n, &m, &tpone, U, &iu, &ju, descU, A, &ia, &ja, descA,
	       		&tzero, VTcpy, &ivtcpy, &jvtcpy, descVTcpy ); 

	for (i=1;i<min_mn+1;i++){
		prow = indxg2p_( &i, &nbA, &izero, &izero, &nprow );
		localcol = indxg2l_( &i, &nbA, &izero, &izero, &nprow );
		invStemp = 1/S[i-1];
		if( myrow==prow )
			sscal_( &nqA, &invStemp, &(VTcpy[localcol-1]), &mpVT );
	}

	verif_repres_VN = verif_orthogonality(min_mn,n,VTcpy, ivtcpy, jvtcpy, descVTcpy);

	free(VTcpy);

	return verif_repres_VN;
}


	int driver_psgesvd( char jobU, char jobVT, int m, int n, float *A, int ia, int ja, int *descA,
		float *S_NN, float *U_NN, int iu, int ju, int *descU, float *VT_NN, int ivt, int jvt, int *descVT,
	       	double *MPIelapsedNN){

	float *Acpy=NULL, *work=NULL;
	int lwork;
/**/
	int ione=1;
/**/
	int nprow, npcol, myrow, mycol;
	int mpA, nqA;
	int ictxt, nbA, rsrcA, csrcA;
	int ctxt_ = 1, nb_ = 5, rsrc_ = 6, csrc_ = 7;
	int infoNN;

	double MPIt1, MPIt2;
/**/
	ictxt = descA[ctxt_];
	Cblacs_gridinfo( ictxt, &nprow, &npcol, &myrow, &mycol );

	nbA = descA[nb_]; rsrcA = descA[rsrc_] ; csrcA = descA[csrc_] ;

	mpA    = numroc_( &m     , &nbA, &myrow, &rsrcA, &nprow );
	nqA    = numroc_( &n     , &nbA, &mycol, &csrcA, &npcol );

	Acpy = (float *)calloc(mpA*nqA,sizeof(float)) ;
	if (Acpy==NULL){ printf("error of memory allocation Acpy on proc %dx%d\n",myrow,mycol); exit(0); }

      	pslacpy_( "All", &m, &n, A, &ione, &ione, descA, Acpy, &ione, &ione, descA );


	work = (float *)calloc(1,sizeof(float)) ;
	if (work==NULL){ printf("error of memory allocation for work on proc %dx%d (1st time)\n",myrow,mycol); exit(0); }

	lwork=-1;

	psgesvd_( &jobU, &jobVT, &m, &n, Acpy, &ione, &ione, descA,
	       	S_NN, U_NN, &ione, &ione, descU, VT_NN, &ione, &ione, descVT,
	       	work, &lwork, &infoNN);

	lwork = (int) (work[0]);
	free(work);

	work = (float *)calloc(lwork,sizeof(float)) ;
	if (work==NULL){ printf("error of memory allocation work on proc %dx%d\n",myrow,mycol); exit(0); }
/**/		
	MPIt1 = MPI_Wtime();
/**/
	psgesvd_( &jobU, &jobVT, &m, &n, Acpy, &ione, &ione, descA,
	       	S_NN, U_NN, &ione, &ione, descU, VT_NN, &ione, &ione, descVT,
	       	work, &lwork, &infoNN);
/**/
	MPIt2 = MPI_Wtime();
	(*MPIelapsedNN)=MPIt2-MPIt1;
/**/
	free(work);
	free(Acpy);
	return infoNN;
}

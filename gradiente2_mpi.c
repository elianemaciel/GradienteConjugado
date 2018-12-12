/*###########################################################################
	AUTOR:	ELIANE MACIEL
			FATIMA
    VERSÃO: 2
    DESENVOLVIMENTO PARALELO
############################################################################*/

#include <sys/time.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "hb_io.h"
#include <math.h>
#include <mpi.h>
#include <time.h>
#define NTHREADS 2

void imprimeResultado(double *resultado, int n, double inicio, double fim){
    /*
    Função que faz a impressão dos vetores resultados
    */
    printf("\nResultado x:\n\n");
	int i;
    for(i=0; i<n; i++){
        printf("%.16f\n",resultado[i]);
    }
    printf("\nTempo = %f\n", fim-inicio);
}

void imprimeProvaReal(double *r, int *rowind, double *values, int *colptr, double *x, int n){
    /*
    Função que faz a impressão de A.x
    */
    int i, coluna;
    for (i = 0; i < n; i++)
        r[i] = 0;
    coluna = -1;
    i = 0;
    while(rowind[i] != NULL){
        if(i + 1 == colptr[coluna + 1])
            coluna++;
        r[rowind[i]-1] += values[i] * x[coluna];
        if ((rowind[i] - 1) != coluna){
                r[coluna] += values[i] * x[rowind[i] - 1];
            }
        i++;
    }
    printf("\nAx = b \n");
    for(i=0;i<n;i++){
        printf("%.16f\n", r[i]);
    }
    printf("\n");
}

void copiaVetor(double *vetor, double *copia, int n){
    /*
    Função que copia valores de um vetor para outro vetor
    */
	int i;
	for (i = 0; i < n; i++){
	    copia[i] = vetor[i];
	}
}

void multiplicacaoVetorValor(double *vetor1, double num, double *res, int n){
    /*
    Função que faz a multiplicação Vetor por um valor
    */
	int i;
	for(i = 0 ; i < n ; i++){
		res[i] = vetor1[i] * num;
	}
}

void somaVetorValor(double *vetor1, double num, double *res, int n){
    /*
    Função que faz a soma de um vetor por um valor
    */
	int i;
	for(i = 0 ; i < n ; i++){
		res[i] = vetor1[i] + num;
	}
}

void somaVetorVetor(double *vetor1, double *vetor2, double *res, int n){
    /*
    Função que faz a soma de dois vetores
    */
	int i;
	for(i = 0 ; i < n ; i++){
		res[i] = vetor1[i] + vetor2[i];
	}
}

void subtracaoVetorVetor(double *vetor1, double *vetor2, double *res, int n){
    /*
    Função que faz a subtração de dois vetores
    */
	int i;
	for(i = 0 ; i < n ; i++){
		res[i] = vetor1[i] - vetor2[i];
	}
}

void geraVetor(double *b, int ncol){
    /*
    Função que faz a geração de valores b
    */
    int i;
	for( i=0 ; i<ncol ; i++ ){
		b[i] = 1;
	}
}

void multiplicacaoMatrizVetor(int id, int n, int np, double *qlocal, double *values, double *d, int *colptr, int *rowind){
    /*
    Função que faz a multiplicação da matriz pelo vetor
    */
    int i,coluna;
    for(coluna = id;coluna<n;coluna+=np){
        i = colptr[coluna]-1;

        do{
            qlocal[rowind[i]-1] += values[i] * d[coluna];
            if ((rowind[i] - 1) != coluna){
                qlocal[coluna] += values[i] * d[rowind[i] - 1];
            }
            i++;
        }while(i < colptr[coluna+1]-1);
    }
}

void multiplicacaoMatrizVetorDois(int id, int n, int np, double *rlocal, double *values, double *x, int *colptr, int *rowind){
    /*
    Função que faz a multiplicação da matriz pelo vetor
    */
    int i, coluna;
    for(coluna = id; coluna<n ;coluna+=np){
        i = colptr[coluna]-1;
        do{
            rlocal[rowind[i]-1] -= values[i] * x[coluna];
            if ((rowind[i] - 1) != coluna){
                rlocal[coluna] -= values[i] * x[rowind[i] - 1];
            }
            i++;
        }while(i < colptr[coluna+1]-1);
    }
}

void gradienteConjugado(double *values, int *colptr, int *rowind, double *b, int ncol, int argc, char *argv[]){
    /*
    Função que faz o gerenciamento do algoritmo gradiente conjugado
    */
    int imax = 1000, coluna, id, np, a = 1, i;
    double erro = 0.00001;
    double *x, *r, *d, *q, *resultado, *qlocal, *rlocal;
    double dq, sigma_novo = 0, sigma0, sigma_velho, alpha, beta, inicio, fim;

    qlocal = (double*)malloc(ncol*sizeof(double));
    rlocal = (double*)malloc(ncol*sizeof(double));

    x = (double*)malloc(ncol*sizeof(double));
    r = (double*)malloc(ncol*sizeof(double));
    d = (double*)malloc(ncol*sizeof(double));
    q = (double*)malloc(ncol*sizeof(double));
    resultado = (double *)malloc(ncol*sizeof(double));

    // x = zeros(n,1);
    for (i = 0; i < ncol; i++)
        x[i] = 0;

    // r = b - A * x;
    copiaVetor(b, r, ncol);

    // d = r;
    copiaVetor(r, d, ncol);

    // sigma_novo = r' * r;
    for(i = 0; i < ncol; i++){
        sigma_novo += r[i] * r[i];
    }

    sigma0 = sigma_novo;

    // MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &np);
    inicio = MPI_Wtime();

    while (a < imax && sigma_novo > erro * erro * sigma0){
        for(i = 0; i < ncol; i++){
            q[i] = 0;
            qlocal[i] = 0;
            rlocal[i] = 0;
        }

        // q = A * d;
		MPI_Barrier(MPI_COMM_WORLD);

        multiplicacaoMatrizVetor(id, ncol, np, qlocal, values, d, colptr, rowind);

		MPI_Allreduce(qlocal,q,ncol,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

		MPI_Barrier(MPI_COMM_WORLD);

        // alpha = sigma_novo/(d' * q);
        dq = 0;
        for(i = 0; i < ncol; i++){
            dq += d[i] * q[i];
        }
        alpha = sigma_novo/dq;

        // x = x + alpha * d;
        multiplicacaoVetorValor( d, alpha, resultado, ncol);
        somaVetorVetor(x, resultado, x, ncol);

        if(a % 50 == 0){
            //r = b - A * x;
            copiaVetor(b, r, ncol);

            MPI_Barrier(MPI_COMM_WORLD);

            multiplicacaoMatrizVetorDois(id, ncol, np, rlocal, values, x, colptr, rowind);

            MPI_Allreduce(rlocal,r,ncol,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

            for(i=0;i<ncol;i++){
	   			r[i] += b[i];
			}
        }
        else{
            // r = r - alpha * q;
            multiplicacaoVetorValor( q, alpha, resultado, ncol);
            subtracaoVetorVetor(r, resultado, r,  ncol);
        }

        sigma_velho = sigma_novo;

        // sigma_novo = r' * r;
        sigma_novo = 0;
        for(i = 0; i < ncol; i++){
            sigma_novo += r[i] * r[i];
        }

        beta = sigma_novo / sigma_velho;

        // d = r + beta * d;
        multiplicacaoVetorValor(d, beta,resultado, ncol);
        somaVetorVetor(r, resultado, d, ncol);
        a++;
    }
    fim = MPI_Wtime();
    MPI_Finalize();
    if ( id == 0 ){
        imprimeResultado(x, ncol, inicio, fim);
        imprimeProvaReal(r, rowind, values, colptr, x, ncol);
    }
}


int main (int argc, char *argv[]) {
    double *b;
    int *colptr = NULL, indcrd;
    char *indfmt = NULL;
    FILE *input, *arq, *file_matriz;
    char *key = NULL;
    char *mxtype = NULL;
    int ncol, neltvl, nnzero, nrhs, nrhsix, nrow;
    int ptrcrd;
    char *ptrfmt = NULL;
    int rhscrd;
    char *rhsfmt = NULL, *rhstyp = NULL, *title = NULL;
    int *rowind = NULL;
    int totcrd, valcrd;
    char *valfmt = NULL;
    double *values = NULL;

    // input = fopen("entradas/matriz/bcsstk11.rsa", "r");
    input = fopen("entradas/matriz/matrizMenor.rsa", "r");

    if ( input == NULL ){
        printf("Erro ao abrir o arquivo\n");
        return 0;
    }

    hb_header_read(input, &title, &key, &totcrd, &ptrcrd, &indcrd,
        &valcrd, &rhscrd, &mxtype, &nrow, &ncol, &nnzero, &neltvl, &ptrfmt,
        &indfmt, &valfmt, &rhsfmt, &rhstyp, &nrhs, &nrhsix
    );

    colptr = ( int * ) malloc ( ( ncol + 1 ) * sizeof ( int ) );

    if ( mxtype[2] == 'A' )
    {
      rowind = ( int * ) malloc ( nnzero * sizeof ( int ) );
      values = ( double * ) malloc ( nnzero * sizeof ( double ) );
    }
    else if ( mxtype[2] == 'E' )
    {
      rowind = ( int * ) malloc ( neltvl * sizeof ( int ) );
      values =  ( double * ) malloc ( neltvl * sizeof ( double ) );
    }
    else
    {
      printf ( "\n" );
      printf ( "TEST05 - Warning!\n" );
      printf ( "  Illegal value of MXTYPE character 3.\n" );
    }


    hb_structure_read ( input, ncol, mxtype, nnzero, neltvl,
      ptrcrd, ptrfmt, indcrd, indfmt, colptr, rowind );


    hb_values_read ( input, valcrd, mxtype, nnzero, neltvl, valfmt, values );

    fclose ( input );

    int i=0;

    b = (double*)malloc(ncol*sizeof(double));

    geraVetor(b, ncol);

    gradienteConjugado(values,colptr,rowind,b,ncol, argc, argv);

    free ( colptr );
    free ( rowind );
    free ( values );


    return 0;

}

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
    printf("\n\nResultado:\n\n");
	int i;
    for(i=0;i<n;i++){
        printf("%f\n",resultado[i]);
    }
    printf("Tempo = %f\n", fim-inicio);
}

void copiaVetor(double *vetor, double *copia, int n){
	int i;
	for (i = 0; i < n; i++){
	    copia[i] = vetor[i];
	}
}

void multiplicacaoVetorValor(double *vetor1, double num, double *res, int n){
	int i;
	for(i = 0 ; i < n ; i++){
		res[i] = vetor1[i] * num;
	}
}

void somaVetorValor(double *vetor1, double num, double *res, int n){
	int i;
	for(i = 0 ; i < n ; i++){
		res[i] = vetor1[i] + num;
	}
}

void somaVetorVetor(double *vetor1, double *vetor2, double *res, int n){
	int i;
	for(i = 0 ; i < n ; i++){
		res[i] = vetor1[i] + vetor2[i];
	}

}

void subtracaoVetorVetor(double *vetor1, double *vetor2, double *res, int n){
	int i;
	for(i = 0 ; i < n ; i++){
		res[i] = vetor1[i] - vetor2[i];
	}
}

void geraVetor(double *b, int ncol){
    int i;
	srand(time(NULL));

	for( i=0 ; i<ncol ; i++ ){
		b[i] = rand()%10;
		printf("b = %f \n", b[i]);
	}
}

void gradienteConjugado(double *values, int *colptr, int *rowind, double *b, int ncol, int argc, char *argv[]){
    int imax = 1000;
    double erro = 0.00001;
    int a = 1, i;
    double *x, *r, *d, *q, *resultado, *qlocal;
    double dq;
    double sigma_novo = 0, sigma0, sigma_velho;
    double alpha, beta, inicio, fim;
    int coluna;
    int id, np;

    qlocal = (double*)malloc(ncol*sizeof(double));

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


	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &np);
    printf("ID = %d NP = %d\n", id, np);
    inicio = MPI_Wtime();

    while (a < imax && sigma_novo > erro * erro * sigma0)
    {
        for(i = 0; i < ncol; i++){
            q[i] = 0;
            qlocal[i] = 0;
        }

        // q = A * d;
        // multiplicacaoMatrizVetor(values, colptr, rowind, d, q);
        int pptr,col;
		MPI_Barrier(MPI_COMM_WORLD);
        // printf("ID = %d\n", id);
		for(col = id;col<ncol;col+=np){
			pptr = colptr[col]-1;

			do{
				qlocal[rowind[pptr]-1] += values[pptr] * d[col];
				pptr++;
			}while(pptr < colptr[col+1]-1);
		}

		MPI_Allreduce(qlocal,q,ncol,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
		// printf("\nq:\n");
		// for(i=0;i<ncol;i++){
		// 	printf("%f\t",q[i]);
		// }
		// printf("\n");

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

            coluna = -1;
            i = 0;
            while(rowind[i] != NULL){
                if(i + 1 == colptr[coluna + 1])
                    coluna++;
                r[rowind[i]-1] -= values[i] * x[coluna];
                i++;
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
        // printf("Tempo = %f\n", fim - inicio);
    }
}


int main (int argc, char *argv[]) {
    double *b;
    int *colptr = NULL;
    int indcrd;
    char *indfmt = NULL;
    FILE *input, *arq, *file_matriz;
    char *key = NULL;
    int khi;
    int klo;
    char *mxtype = NULL;
    int ncol;
    int neltvl;
    int nnzero;
    int nrhs;
    int nrhsix;
    int nrow;
    int ptrcrd;
    char *ptrfmt = NULL;
    int rhscrd;
    char *rhsfmt = NULL;
    char *rhstyp = NULL;
    int *rowind = NULL;
    char *title = NULL;
    int totcrd;
    int valcrd;
    char *valfmt = NULL;
    double *values = NULL;

    int id, nproc, resto = 0, ult_linha;

    input = fopen("entradas/matriz/bcsstk05.rsa", "r");
    // input = fopen("entradas/matriz/matrizMenor.rsa", "r");

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

    // arq = fopen("entradas/vetor/vetor.txt", "r");
    // arq = fopen("entradas/vetor/vetorMenor.txt", "r");

    // i = 0;
    // char linha[3];
    // char *result;
    // while (!feof(arq))
    // {
    // // Lê uma linha (inclusive com o '\n')
    //     result = fgets(linha, 3, arq);  // o 'fgets' lê até 3 caracteres ou até o '\n'
    //     if (result){ // Se foi possível ler
    //         if(linha != NULL){
    //             b[i] = atof(linha);
    //         }
    //     }

    //     i++;
    // }
    // fclose(arq);

    gradienteConjugado(values,colptr,rowind,b,ncol, argc, argv);

    free ( colptr );
    free ( rowind );
    free ( values );


    return 0;

}
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
#include <omp.h>
#include <time.h>


#define NTHREADS 2

void imprimeResultado(double *resultado, int n){
    printf("\n\nResultado:\n\n");
	int i;
    for(i=0;i<n;i++){
        printf("%f\n",resultado[i]);
    }
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
    int imax = 2;
    double erro = 0.00001;
    int a = 1, i;
    double *x, *r, *d, *q, *resultado, *qlocal;
    double dq;
    double sigma_novo = 0, sigma0, sigma_velho;
    double alpha, beta;
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


	// MPI_Init(&argc, &argv);
	// MPI_Comm_rank(MPI_COMM_WORLD, &id);
	// MPI_Comm_size(MPI_COMM_WORLD, &np);

    while (a < imax && sigma_novo > erro * erro * sigma0)
    {
        for(i = 0; i < ncol; i++){
            q[i] = 0;
            qlocal[i] = 0;
        }

        // q = A * d;

        int id, i, ini, fim,parte;

        int j = 0;
        int pptr,col;
        // coluna e j não são private porque compartilham nas threds
        #pragma omp parallel private (id,ini,fim,parte, pptr)
        {
            parte = 9 / omp_get_num_threads(); // Total de elementos dividido pelas threds
            id = omp_get_thread_num( );
            ini = parte*id;
            fim =ini+parte;

            printf("ID = %d, ini =%d fim = %d\n", id, ini, fim);
            for(col=ini;col<fim;col++){
    			pptr = colptr[col]-1;
                printf("pptr = %d col = %d \n", pptr, col);

    			do{
    				qlocal[rowind[pptr]-1] += values[pptr] * d[col];
                    printf("values[%d] = %f d[%d] = %f \n", pptr, values[pptr], col, d[col]);
    				pptr++;
    			}while(pptr < colptr[col+1]-1);
    		}
            #pragma omp critical
            {
                somaVetorVetor(qlocal, q, q, ncol);
            }
        }

        printf("\nq:\n");
		for(i=0;i<ncol;i++){
			printf("%f\t",q[i]);
		}
		printf("\n");

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
    // MPI_Finalize();
    imprimeResultado(x, ncol);
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

    omp_set_num_threads(NTHREADS);
    // input = fopen("entradas/matriz/bcsstruc2.data", "r");
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

    // geraVetor(b, ncol);

    // arq = fopen("entradas/vetor/vetor.txt", "r");
    arq = fopen("entradas/vetor/vetorMenor.txt", "r");

    i = 0;
    char linha[3];
    char *result;
    while (!feof(arq))
    {
    // Lê uma linha (inclusive com o '\n')
        result = fgets(linha, 3, arq);  // o 'fgets' lê até 3 caracteres ou até o '\n'
        if (result){ // Se foi possível ler
            if(linha != NULL){
                b[i] = atof(linha);
            }
        }

        i++;
    }
    fclose(arq);

    gradienteConjugado(values,colptr,rowind,b,ncol, argc, argv);

    free ( colptr );
    free ( rowind );
    free ( values );


    return 0;

}

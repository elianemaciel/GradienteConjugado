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

void multiplicacaoMatrizVetor(double *values, int *colptr, int *rowind, double *vetor, double *resultado){
    int coluna, i;
    coluna = -1;
    i = 0;
    while(rowind[i] != NULL){
        if(i + 1 == colptr[coluna + 1]){
            coluna++;
        }
        resultado[rowind[i] - 1] += values[i] * vetor[coluna];
        i++;
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
    // Inicializa vetor com números randomicos
	srand(time(NULL));
    int i;
	for( i=0 ; i<ncol ; i++ ){
		b[i] = rand()%10;
		printf("b = %f \n", b[i]);
	}
}

// void multiplicacao_vetor_matriz(int *rowind, int *colptr, double *values, double *q, double *d, int N){
//     //q = A * d;
//     // Ainda não esta correto o paralelo Cada execução um valor diferente
//     // Verificar como fazer a parte espelhada
//
//     double *v;
//
//     v = (double*)malloc(9*sizeof(double));
//
//     for(int j = 0;j < 9;j++){
//         v[j] = 0;
//     }
//     int id, i, ini, fim,parte, coluna;
//     double prodparc;
//
//     coluna = -1;
//     int j = 0;
//     // coluna e j não são private porque compartilham nas threds
//     #pragma omp parallel private (id,ini,fim,parte)
//     {
//         parte = 9 / omp_get_num_threads(); // Total de elementos dividido pelas threds
//         id = omp_get_thread_num( );
//         ini = parte*id;
//         fim =ini+parte;
//
//         printf("ID = %d\n", id);
//         for (i=ini;i<fim;i++){
//             if(j + 1 == colptr[coluna + 1]){
//                 coluna++;
//             }
//
//             // Verificar como colocar no vetor
//             v[i] += values[j] * d[coluna];
//             printf("values[%d] = %f d[%d] = %f v[%d] = %f \n", j, values[j], coluna, d[coluna], rowind[j]-1,q[rowind[j] - 1]);
//             j++;
//         }
//         #pragma omp critical
//         {
//             // coloca a soma no vetor principal de resultados
//             // q[rowind[i] - 1] += values[i] * d[coluna];
//             // prod = prod + prodparc;
//             for (i=ini;i<fim;i++){
//                 q[rowind[i] - 1] += v[i];
//                 printf("q[%d] = %f \n", i, q[i]);
//             }
//         }
//     }
// }


void gradienteConjugado(double *values, int *colptr, int *rowind, double *b, int ncol){
    int imax = 1000;
    double erro = 0.00001;
    int a = 1, i;
    double *x, *r, *d, *q, *resultado;
    double dq;
    double sigma_novo = 0, sigma0, sigma_velho;
    double alpha, beta;
    int coluna;

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

    while (a < imax && sigma_novo > erro * erro * sigma0)
    {
        for(i = 0; i < ncol; i++){
            q[i] = 0;
        }

        // q = A * d;
        multiplicacaoMatrizVetor(values, colptr, rowind, d, q);

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

    // omp_set_num_threads(NTHREADS);


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

    gradienteConjugado(values,colptr,rowind,b,ncol);

    free ( colptr );
    free ( rowind );
    free ( values );


    return 0;

}

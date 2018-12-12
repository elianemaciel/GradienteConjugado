/*###########################################################################
	AUTOR:	ELIANE MACIEL
			FATIMA
    VERSÃO: 1
############################################################################*/

#include <sys/time.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "hb_io.h"
#include <math.h>
#include <time.h>

void imprimeResultado(double *resultado, int n, clock_t inicio, clock_t fim){
    printf("\n\nResultado:\n\n");
    int i;
    clock_t t;
    for(i=0;i<n;i++){
        printf("%.20f\n",resultado[i]);
    }
    t = fim - inicio;
    double time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds

    printf("\nTempo = %f\n", time_taken);
}

void imprimeProvaReal(double *r, int *rowind, double *values, int *colptr, double *x, int n){
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

void gradienteConjugado(double *values, int *colptr, int *rowind, double *b, int n){
    int imax = 1000;
    double erro = 0.0000001;
    int a = 1, i;
    double *x, *r, *d, *q;
    double dq;
    double sigma_novo = 0, sigma0, sigma_velho;
    double alpha, beta;
    int coluna;

    x = (double*)malloc(n*sizeof(double));
    r = (double*)malloc(n*sizeof(double));
    d = (double*)malloc(n*sizeof(double));
    q = (double*)malloc(n*sizeof(double));

    // x = zeros(n,1);
    for (i = 0; i < n; i++)
        x[i] = 0;

    // r = b - A * x;
    for(i = 0; i < n; i++){
        r[i] = b[i];
    }

    // d = r;
    for(i = 0;i < n;i++){
        d[i] = r[i];
    }

    // sigma_novo = r' * r;
    for(i = 0; i < n; i++){
        sigma_novo += r[i] * r[i];
    }

    sigma0 = sigma_novo;

    clock_t begin=clock();

    while (a < imax && sigma_novo > erro * erro * sigma0)
    {
        for(i = 0; i < n; i++){
            q[i] = 0;
        }

        //q = A * d;
        coluna = -1;
        i = 0;
        while(rowind[i] != NULL){
            if(i + 1 == colptr[coluna + 1]){
                coluna++;
            }
            q[rowind[i] - 1] += values[i] * d[coluna];
            if ((rowind[i] - 1) != coluna){
                q[coluna] += values[i] * d[rowind[i] - 1];
            }
            i++;
        }

        // alpha = sigma_novo/(d' * q);
        dq = 0;
        for(i = 0; i < n; i++){
            dq += d[i] * q[i];
        }
        alpha = sigma_novo/dq;

        // x = x + alpha * d;
        for(i = 0; i < n; i++){
            x[i] += alpha * d[i];
        }

        if(a % 50 == 0){
            //r = b - A * x;
            for(i = 0; i < n; i++){
                r[i] = b[i];
            }

            coluna = -1;
            i = 0;
            while(rowind[i] != NULL){
                if(i + 1 == colptr[coluna + 1])
                    coluna++;
                r[rowind[i]-1] -= values[i] * x[coluna];
                if ((rowind[i] - 1) != coluna){
                    r[coluna] -= values[i] * x[rowind[i] - 1];
                }
                i++;
            }
        }
        else{
            // r = r - alpha * q;
            for(i = 0; i < n; i++){
                r[i] += - alpha * q[i];
            }
        }

        sigma_velho = sigma_novo;

        // sigma_novo = r' * r;
        sigma_novo = 0;
        for(i = 0; i < n; i++){
            sigma_novo += r[i] * r[i];
        }


        beta = sigma_novo / sigma_velho;

        // d = r + beta * d;
        for(i = 0; i < n; i++){
            d[i] = r[i] + beta * d[i];
        }
        a++;
    }
    clock_t end=clock();
    imprimeResultado(x, n, begin, end);
    imprimeProvaReal(r, rowind, values, colptr, x, n);
}

void geraVetor(double *b, int ncol){
    int i;
    srand(time(NULL));

    for( i=0 ; i<ncol ; i++ ){
        b[i] = 1;
    }
}


int main (int argc, char *argv[]) {
    double *b;
    int *colptr = NULL;
    int indcrd;
    char *indfmt = NULL;
    FILE *input, *arq, *file_matriz;
    char *key = NULL;
    int khi, klo;
    char *mxtype = NULL;
    int ncol, neltvl, nnzero, nrhs, nrhsix, nrow, ptrcrd;
    char *ptrfmt = NULL;
    int rhscrd;
    char *rhsfmt = NULL, *rhstyp = NULL, *title = NULL, *valfmt = NULL;
    int *rowind = NULL;
    int totcrd, valcrd;
    double *values = NULL;

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

    // Vetor colptr
    b = (double*)malloc(ncol*sizeof(double));
    geraVetor(b, ncol);

    gradienteConjugado(values,colptr,rowind,b,ncol);

    free ( colptr );
    free ( rowind );
    free ( values );

    return 0;

}

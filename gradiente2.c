#include <sys/time.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "hb_io.h"
#include <math.h>
#include<omp.h>

#define NTHREADS 2

void imprimeResultado(double *resultado, int n){
    printf("\n\nResultado:\n\n");
    for(int i=0;i<n;i++){
        printf("%f\n",resultado[i]);
    }
}

void multiplicacao_vetor_matriz(int *rowind, int *colptr, double *values, double *q, double *d, int N){
    //q = A * d;
    // Ainda não esta correto o paralelo Cada execução um valor diferente
    // Verificar como fazer a parte espelhada

    double *v;

    v = (double*)malloc(9*sizeof(double));

    for(int j = 0;j < 9;j++){
        v[j] = 0;
    }
    int id, i, ini, fim,parte, coluna;
    double prodparc;
 
    coluna = -1;
    int j = 0;
    // coluna e j não são private porque compartilham nas threds 
    #pragma omp parallel private (id,ini,fim,parte)
    {
        parte = 9 / omp_get_num_threads(); // Total de elementos dividido pelas threds
        id = omp_get_thread_num( );
        ini = parte*id;
        fim =ini+parte;

        printf("ID = %d\n", id);
        for (i=ini;i<fim;i++){
            if(j + 1 == colptr[coluna + 1]){
                coluna++;
            }
             
            // Verificar como colocar no vetor
            v[i] += values[j] * d[coluna];
            printf("values[%d] = %f d[%d] = %f v[%d] = %f \n", j, values[j], coluna, d[coluna], rowind[j]-1,q[rowind[j] - 1]);
            j++;
        }
        #pragma omp critical
        {
            // coloca a soma no vetor principal de resultados
            // q[rowind[i] - 1] += values[i] * d[coluna];
            // prod = prod + prodparc;
            for (i=ini;i<fim;i++){
                q[rowind[i] - 1] += v[i];
                printf("q[%d] = %f \n", i, q[i]);
            }
        }
    }
}

void gradienteConjugado(double *values, int *colptr, int *rowind, double *b, int n){
    int imax = 1000;
    double erro = 0.00001;
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

    // multiplicacao_vetor_matriz(rowind,colptr, values, q, d, n);

    while (a < imax && sigma_novo > erro * erro * sigma0)
    {
        for(i = 0; i < n; i++){
            q[i] = 0;
        }

        multiplicacao_vetor_matriz(rowind,colptr, values, q, d, n);

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
    printf("\n\nResultado:\n\n");
    for(int i=0;i<n;i++){
        printf("%.2f\n",x[i]);
    }
    // imprimeResultado(x, n);
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

    // arq = fopen("entradas/vetor/vetor.txt", "r");
    arq = fopen("entradas/vetor/vetorMenor.txt", "r");

    i = 0;
    char linha[3];
    char *result;
    while (!feof(arq))
    {
    // Lê uma linha (inclusive com o '\n')
        result = fgets(linha, 3, arq);  // o 'fgets' lê até 3 caracteres ou até o '\n'
        //printf("Linha lida %s",result);
        if (result){ // Se foi possível ler
            printf("Linha %d : %s",i,linha);
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

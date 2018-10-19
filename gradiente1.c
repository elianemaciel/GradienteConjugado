#include <sys/time.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "hb_io.h"
#include <math.h>

void imprimeResultado(double *resultado, int n){
    printf("\n\nResultado:\n\n");
    for(int i=0;i<n;i++){
        printf("%f\n",resultado[i]);
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
    imprimeResultado(x, n);
}


int main (int argc, char *argv[]) {
    double *b;
    int *colptr = NULL;
    int indcrd;
    char *indfmt = NULL;
    FILE *input;
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

    if (argc < 2){
        fprintf(stderr, "%s < arquivo matriz >\n", argv[0]);
        return 0;
    }

    input = fopen(argv[1], "r");

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

    printf ( "\n" );
    printf ( "  '%s'\n", title );
    printf ( "  KEY =    '%s'\n", key );
    printf ( "\n" );
    printf ( "  NROW =   %d\n", nrow );
    printf ( "  NCOL =   %d\n", ncol );
    printf ( "  NNZERO = %d\n", nnzero );
    printf ( "  NELTVL = %d\n", neltvl );

    int i=0;

    // Vetor colptr
    for(i=0;i<nrow;i++){
        printf ( "  colptr =   %d\n", colptr[i] );
    }

    printf("\n");

    // Matriz rowind
    for(i=0; i<nnzero;i++){
        printf ( "  rowind =   %d\n", rowind[i] );
    }

    // Valores
    printf("\n");
    for(i=0; i<nnzero;i++){
        printf ( "  values =   %.2f\n", values[i] );
    }

    b = (double*)malloc(ncol*sizeof(double));
    printf("\nInforme o Vetor:\n");
    for(i=0;i<ncol;i++){
        scanf("%lf", &b[i]);
    }               

    gradienteConjugado(values,colptr,rowind,b,ncol);

    free ( colptr );
    free ( rowind );
    free ( values );
    
    return 0;

}

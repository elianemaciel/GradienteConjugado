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
    int c,i,j,k;
    double *x,*r,*d,sum;
    double sigma_novo,sigma0,sigma_velho;
    double *q;
    double *aux1;
    double alpha,beta;
    
    x = (double*)malloc(n*sizeof(double)); 
    r = (double*)malloc(n*sizeof(double));    
    d = (double*)malloc(n*sizeof(double)); 
    q = (double*)malloc(n*sizeof(double));   

    for (int i=0;i<n;i++)
        x[i] = 0;

    for(i=0;i<n;i++){
        r[i] = b[i];
    }
    
    for(i=0;i<n;i++){
        d[i] = r[i];
    }
    
    sigma_novo = 0;
    for(i=0;i<n;i++){
        sigma_novo += r[i]*r[i]; 
    }
    sigma0 = sigma_novo;    
    c=1;
    while (c < imax && sigma_novo > erro * erro * sigma0)
    {
        for(i=0;i<n;i++)
            q[i]=0;
            
        int col=-1;
        i = 0;
        while(rowind[i]!=NULL){
            if(i+1==colptr[col+1]){ 
                col++;
            }
            q[rowind[i]-1] += values[i] * d[col];   
            i++;        
        }

        sum = 0;
        for(i=0;i<n;i++){
            sum += d[i]*q[i];
        }
        alpha = sigma_novo/sum;
        
        for(i=0;i<n;i++){
            x[i] += alpha *d[i];
        }   
        
        if(c%50 == 0){                    
            for(i=0;i<n;i++)
                r[i] = b[i];
            
            int col=-1;
            i = 0;
            while(rowind[i]!=NULL){
                if(i+1==colptr[col+1])
                    col++;
                r[rowind[i]-1] -= values[i] * x[col];   
                i++;        
            }
            
        }
        else{
            for(i=0;i<n;i++){
                r[i] += - alpha * q[i];
            }
        }
        
        sigma_velho = sigma_novo;
        sigma_novo = 0;
        for(i=0;i<n;i++){
            sigma_novo += r[i]*r[i]; 
        }
        
        beta = sigma_novo / sigma_velho;
        for(i=0;i<n;i++){
            d[i] = r[i] + beta * d[i]; 
        }
        c++;
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
    for(i=0;i<ncol;i++){
        b[i] = colptr[i];
    }                

    gradienteConjugado(values,colptr,rowind,b,ncol);

    free ( colptr );
    free ( rowind );
    free ( values );
    
    return 0;

}

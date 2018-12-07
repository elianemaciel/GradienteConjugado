#include <sys/time.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "hb_io.h"
#include <math.h>
#include <mpi.h>

void imprimeResultado(double *resultado, int n){
    printf("\n\nResultado:\n\n");
    for(int i=0;i<n;i++){
        printf("%f\n",resultado[i]);
    }
}

void multiplicacao_vetor_matriz(double *values, int *colptr, int *rowind, double *q, double *d, int i, int coluna, ){
    int *displs, *scounts, *sdiv, *sdivC;
    displs = (int *)malloc(np*sizeof(int)); 
    scounts = (int *)malloc(np*sizeof(int)); 
    /**************************************************************************
        cria os vetores com as informações
        displs - vetor com o valor de onde começa cada processo dentro do col_ptr
        scounts - vetor com a quantidade que será enviada para cada processo
    **************************************************************************/
    for (i=0; i<np; ++i) { 
        displs[i] = i*(div); 
        scounts[i] = div+1; 
    } 

    if(id == 0)
        printf("enviando col_ptr\n");   
    //envia o vetor dividido utilizando um offset
    MPI_Scatterv(colptr,scounts, displs,MPI_INT,colptr,div+1,MPI_INT,0,MPI_COMM_WORLD);


    /***********************************************************************************
    Para o servidor:
        sdiv - tem o valor do col_ptr que inicia  o processo i 
        sdivC - tem a quantidade que vai ser enviada para os clientes, porem o cálculo 
                da ultima é feito de forma diferente pois pode ter menos o col_ptr

    Para Cliente:
        sdiv - tem o valor do col_ptr que inicia o processo, sempre busca do 0, pois como foi
                dividido cada processo vai sempre iniciar no col_ptr[0]
        sdivC - tem a quantidade que vai ser recebida, porem o cálulo da última também vai
                ser diferente pois vai ter o delimitador o E.
    ***********************************************************************************/
    if(id == 0 ){       
        sdiv = (int *)malloc(np*sizeof(int)); 
        sdivC = (int *)malloc(np*sizeof(int)); 
        for (i=0; i<np; ++i) { 
            sdiv[i] = colptr[i*(div)];
            if(i != np-1)
                sdivC[i] = colptr[(i*(div))+(div)] - sdiv[i];
            else
                sdivC[i] = colptr[n] - sdiv[i];
        } 
    }
    else{       
        sdiv = (int *)malloc(1*sizeof(int)); 
        sdivC = (int *)malloc(1*sizeof(int));
        sdiv[0] = colptr[0];
        if(id == np-1)
            sdivC[0] =  E - colptr[0];
        else
            sdivC[0] = colptr[div] - colptr[0];
            
        row_ind = (int *) malloc(sizeof(int) * sdivC[0]);
        val = (double *) malloc(sizeof(double) * sdivC[0]);     
        int ini = colptr[0];
        //alterado o valor do col_ptr, pois no row_ind sempre vai começar no 0 para cada processo
        for(i=0; i<=div;i++)
            colptr[i] = colptr[i] - ini;
    }
    
    //Envia o row_ind
    if(id == 0){
        printf("enviando row_ind\n");
        for(i=1;i<np;i++)
            MPI_Send(&row_ind[sdiv[i]],sdivC[i],MPI_INT, i,100,MPI_COMM_WORLD);
    }
    else{
        MPI_Recv(row_ind,sdivC[0],MPI_INT, 0,100,MPI_COMM_WORLD,&s);
    }

    //envia o indicador de limite de linha
    int row_limit = 0;
    if(id == 0){
        row_limit = displs[1];
        for(i=1;i<np-1;i++)
            MPI_Send(&displs[i+1],1,MPI_INT, i,300,MPI_COMM_WORLD);
        displs[np-1]++;
        MPI_Send(&displs[np-1],1,MPI_INT, i,300,MPI_COMM_WORLD);
    }
    else{
        MPI_Recv(&row_limit,1,MPI_INT, 0,300,MPI_COMM_WORLD,&s);
    }
    

    //Envia o vetor val
    if(id == 0){
        printf("enviando valores\n");
        for(i=1;i<np;i++)
            MPI_Send(&val[sdiv[i]],sdivC[i],MPI_DOUBLE, i,200,MPI_COMM_WORLD);
    }
    else{
        MPI_Recv(val,sdivC[0],MPI_DOUBLE, 0,200,MPI_COMM_WORLD,&s);
    }



    // while(rowind[i] != NULL){

    //     if(i + 1 == colptr[coluna + 1]){ 
    //         coluna++;
    //     }
    //     q[rowind[i] - 1] += values[i] * d[coluna];   
    //     i++;        
    // }


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
       
        multiplicacao_vetor_matriz(values,colptr,rowind,q,d,i, coluna);

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

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);


    // if (argc < 2){
    //     fprintf(stderr, "%s < arquivo matriz >\n", argv[0]);
    //     return 0;
    // }

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

    /*printf ( "\n" );
    printf ( "  '%s'\n", title );
    printf ( "  KEY =    '%s'\n", key );
    printf ( "\n" );
    printf ( "  NROW =   %d\n", nrow );
    printf ( "  NCOL =   %d\n", ncol );
    printf ( "  NNZERO = %d\n", nnzero );
    printf ( "  NELTVL = %d\n", neltvl );*/

    int i=0;

    /*
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

    printf("%d\n", ncol);
    */
    b = (double*)malloc(ncol*sizeof(double));

    /*printf("\nInforme o Vetor:\n");
    for(i=0;i<ncol;i++){
        scanf("%lf", &b[i]);
    } */

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

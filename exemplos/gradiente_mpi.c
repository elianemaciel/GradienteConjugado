#include <sys/time.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "hb_io.h"
#include <math.h>
#include <mpi.h>
#include <time.h>

void imprimeResultado(double *resultado, int n){
    printf("\n\nResultado:\n\n");
	int i;
    for(i=0;i<n;i++){
        printf("%f\n",resultado[i]);
    }
}

void multiplicacaoMatrizVetor (double *values, int *colptr, int *rowind, double *vetor, double *resultado, int n, int nc, int id){
	int i,j,l;
	for(i=0;i<id;i++)
		resultado[i] = 0;

	for(i=n,l=0;i<nc;i++,l++){
		for(j=colptr[l];j<colptr[l+1];j++){
			resultado[i] += values[j] * vetor[rowind[j]];
			if( i  != rowind[j]){
				resultado[rowind[j]] += values[j] * vetor[i];
			}
		}
	}
}

float multiplicacaoVetorVetor(double *vetor1, double *vetor2, int n){
	int i;
	float result = 0;
	for (i = 0; i < n; i++){
	    result+= vetor1[i] * vetor2[i];
	}
	return result;
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


int main (int argc, char *argv[]) {
	FILE *input;
	MPI_Status status;
    double *b, *values = NULL, *x, *r, *d, *q, *resultado1 = NULL, *resultado = NULL;
    int ncol, indcrd, neltvl, nnzero, nrhs, nrhsix, nrow, ptrcrd, rhscrd, totcrd, valcrd, id, nProcessos, a, i, acabou, div;
    char *rhsfmt = NULL, *ptrfmt = NULL, *rhstyp = NULL, *key = NULL, *mxtype = NULL, *indfmt = NULL, *title = NULL, *valfmt = NULL;
    int *rowind = NULL, *colptr = NULL, *displs, *scounts, *sdiv, *sdivC;
    double erro = 0.00001, sigmaNovo = 0, sigma0, sigmaVelho, alpha, beta, imax = 1000;;

    input = fopen("entradas/matriz/matrizMenor.rsa", "r");

    if ( input == NULL ){
        printf("Erro ao abrir o arquivo\n");
        return 0;
    }

    MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &nProcessos);

    hb_header_read(input, &title, &key, &totcrd, &ptrcrd, &indcrd,
        &valcrd, &rhscrd, &mxtype, &nrow, &ncol, &nnzero, &neltvl, &ptrfmt,
        &indfmt, &valfmt, &rhsfmt, &rhstyp, &nrhs, &nrhsix
    );

	//envia para todos processos os tamanhos da matriz e vetor
	MPI_Bcast(&nnzero, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&ncol, 1, MPI_INT, 0, MPI_COMM_WORLD);

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

    hb_structure_read ( input, ncol, mxtype, nnzero, neltvl, ptrcrd, ptrfmt, indcrd, indfmt, colptr, rowind );

    hb_values_read ( input, valcrd, mxtype, nnzero, neltvl, valfmt, values);

    fclose ( input );

    b = (double*)malloc(ncol*sizeof(double));
	if(id == 0){
		srand(time(NULL));
		for( i=0 ; i<ncol ; i++ ){
			b[i] = rand()%10;
			printf("b = %f \n", b[i]);
		}
	}

	//divide o tamanho do vetor pelos processos
	div = (1+ncol)/nProcessos;

    x = (double*)malloc(ncol*sizeof(double));
    r = (double*)malloc(ncol*sizeof(double));
    d = (double*)malloc(ncol*sizeof(double));
    q = (double*)malloc(ncol*sizeof(double));
	resultado = (double *)malloc(ncol*sizeof(double));
	displs = (int *)malloc(nProcessos*sizeof(int));
	scounts = (int *)malloc(nProcessos*sizeof(int));

	/**************************************************************************
		cria os vetores com as informações
		displs - vetor com o valor de onde começa cada processo dentro do colptr
		scounts - vetor com a quantidade que será enviada para cada processo
	**************************************************************************/
	for (i=0; i<nProcessos; ++i) {
		displs[i] = i*(div);
		scounts[i] = div+1;
	}

	if(id == 0){
		//envia o vetor dividido utilizando um offset
		MPI_Scatterv(colptr,scounts, displs,MPI_INT,colptr,div+1,MPI_INT,0,MPI_COMM_WORLD);
	}

	/***********************************************************************************
	Para o servidor:
		sdiv - tem o valor do colptr que inicia  o processo i
		sdivC - tem a quantidade que vai ser enviada para os clientes, porem o cálculo
				da ultima é feito de forma diferente pois pode ter menos o colptr

	Para Cliente:
		sdiv - tem o valor do colptr que inicia o processo, sempre busca do 0, pois como foi
				dividido cada processo vai sempre iniciar no colptr[0]
		sdivC - tem a quantidade que vai ser recebida, porem o cálulo da última também vai
				ser diferente pois vai ter o delimitador o nnzero.
	***********************************************************************************/
	if(id == 0 ){
		sdiv = (int *)malloc(nProcessos*sizeof(int));
		sdivC = (int *)malloc(nProcessos*sizeof(int));
		for (i=0; i<nProcessos; ++i) {
			sdiv[i] = colptr[i*(div)];
			if(i != nProcessos-1)
				sdivC[i] = colptr[(i*(div))+(div)] - sdiv[i];
			else
				sdivC[i] = colptr[ncol] - sdiv[i];
		}
	}
	else{
		sdiv = (int *)malloc(1*sizeof(int));
		sdivC = (int *)malloc(1*sizeof(int));
		sdiv[0] = colptr[0];
		if(id == nProcessos-1)
			sdivC[0] =  nnzero - colptr[0];
		else
			sdivC[0] = colptr[div] - colptr[0];

		rowind = (int *) malloc(sizeof(int) * sdivC[0]);
		values = (double *) malloc(sizeof(double) * sdivC[0]);
		int ini = colptr[0];
		//alterado o valor do colptr, pois no rowind sempre vai começar no 0 para cada processo
		for(i=0; i<=div;i++)
			colptr[i] = colptr[i] - ini;
	}

	//Envia o rowind
	if(id == 0){
		for(i=1;i<nProcessos;i++){
			printf("send 1");
			MPI_Send(&rowind[sdiv[i]],sdivC[i],MPI_INT, i,100,MPI_COMM_WORLD);
		}
	}
	else{
		MPI_Recv(rowind,sdivC[0],MPI_INT, 0,100,MPI_COMM_WORLD,&status);
	}

	//envia o indicador de limite de linha
	int row_limit = 0;
	if(id == 0){
		row_limit = displs[1];
		for(i=1;i<nProcessos-1;i++){
			printf("send 2");
			MPI_Send(&displs[i+1],1,MPI_INT, i,300,MPI_COMM_WORLD);
		}
		displs[nProcessos-1]++;
		printf("send 3");
		MPI_Send(&displs[nProcessos-1],1,MPI_INT, i,300,MPI_COMM_WORLD);
	}
	else{
		MPI_Recv(&row_limit,1,MPI_INT, 0,300,MPI_COMM_WORLD,&status);
	}


	//Envia o vetor values
	if(id == 0){
		printf("enviando valores\n");
		for(i=1;i<nProcessos;i++){
			printf("send 4");
			MPI_Send(&values[sdiv[i]],sdivC[i],MPI_DOUBLE, i,200,MPI_COMM_WORLD);
		}
	}
	else{
		MPI_Recv(values,sdivC[0],MPI_DOUBLE, 0,200,MPI_COMM_WORLD,&status);
	}


	// x = zeros(ncol,1);
	x = (double *) malloc(sizeof(double ) * ncol);
	for (i = 0; i < ncol; i++) {
		x[i] = 0;
    	}
	d = (double *) malloc(sizeof(double ) * ncol);


	//r = b - A * x;
	multiplicacaoMatrizVetor(values, colptr, rowind, x, resultado, displs[id], row_limit,ncol);

	if(id == 0){
		for(i=1;i<nProcessos;i++){
			MPI_Recv(resultado1,ncol,MPI_DOUBLE, i,10,MPI_COMM_WORLD,&status);
			somaVetorVetor(resultado,resultado1,resultado,ncol);
		}
		subtracaoVetorVetor(b, resultado, r, ncol);

		//d=r
		copiaVetor(r,d,ncol);

		//sigmaNovo = r' * r;
		sigmaNovo = multiplicacaoVetorVetor(r, r,ncol);
		sigma0 = sigmaNovo;

		acabou = 1;
		a = 1;
		while(a < imax && sigmaNovo > ((erro*erro)*sigma0)){
			for(i=1;i<nProcessos;i++){
				printf("send 5");
				MPI_Send(&acabou,1,MPI_INT, i,20,MPI_COMM_WORLD);
				printf("send 6");
				MPI_Send(d,ncol,MPI_DOUBLE, i,200,MPI_COMM_WORLD);
			}

			// q = A * d;
			multiplicacaoMatrizVetor(values, colptr, rowind, d, q, displs[id], row_limit,ncol);
			for(i=1;i<nProcessos;i++){
				MPI_Recv(resultado1,ncol,MPI_DOUBLE, i,10,MPI_COMM_WORLD,&status);
				somaVetorVetor(q,resultado1,q,ncol);
			}
			//alpha = sigmaNovo/(d' * q);
			alpha = sigmaNovo/multiplicacaoVetorVetor(d, q,ncol);

			//x = x + alpha * d;
			multiplicacaoVetorValor( d, alpha, resultado, ncol);
			somaVetorVetor(x, resultado, x, ncol);

			if(a%50 == 0){
				for(i=1;i<nProcessos;i++){
					printf("send 7");
					MPI_Send(&acabou,1,MPI_INT, i,20,MPI_COMM_WORLD);
					printf("send 8");
					MPI_Send(x,ncol,MPI_DOUBLE, i,200,MPI_COMM_WORLD);
				}
				//r = b - A * x;
				multiplicacaoMatrizVetor(values, colptr, rowind, x, resultado, displs[id], row_limit,ncol);
				for(i=1;i<nProcessos;i++){
					MPI_Recv(resultado1,ncol,MPI_DOUBLE, i,10,MPI_COMM_WORLD,&status);
					somaVetorVetor(resultado,resultado1,resultado,ncol);
				}
				subtracaoVetorVetor(b, resultado, r, ncol);
			}
			else{
				//r = r - alpha * q;
				multiplicacaoVetorValor( q, alpha, resultado, ncol);
				subtracaoVetorVetor(r, resultado, r,  ncol);
			}

			sigmaVelho = sigmaNovo;
			sigmaNovo = multiplicacaoVetorVetor(r, r,ncol);
			double beta = sigmaNovo / sigmaVelho;

			// d = r + beta * d;
			multiplicacaoVetorValor(d, beta,resultado, ncol);
			somaVetorVetor(r, resultado, d, ncol);

			a++;
		}
		// envio de 1 para informar aos processos que terminou.
		acabou=0;
		for(i=1;i<nProcessos;i++){
			printf("send 9");
			MPI_Send(&acabou,1,MPI_INT,i,20,MPI_COMM_WORLD);
		}

		printf("Finalizou com o i = %d\n\n", a);
		imprimeResultado(x, ncol);
	}
	else{
		acabou = 1;
		printf("send 10");
		MPI_Send(resultado,ncol,MPI_DOUBLE, 0,10,MPI_COMM_WORLD);
		MPI_Recv(&ncol,1,MPI_INT, 0,20,MPI_COMM_WORLD,&status);
		while(ncol){
			MPI_Recv(x,ncol,MPI_DOUBLE, 0,200,MPI_COMM_WORLD,&status);
			multiplicacaoMatrizVetor(values, colptr, rowind, x, resultado, displs[id], row_limit,ncol);
			printf("send 11");
			MPI_Send(resultado,ncol,MPI_DOUBLE, 0,10,MPI_COMM_WORLD);
			MPI_Recv(&ncol,1,MPI_INT, 0,20,MPI_COMM_WORLD,&status);
		}
	}

    free ( colptr );
    free ( rowind );
    free ( values );
	free(r);
	free(b);
	free(d);
	free(q);
	free(x);
	free(resultado);
	free(resultado1);

	MPI_Finalize();

    return 0;

}

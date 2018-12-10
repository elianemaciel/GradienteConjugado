#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define N 10

int main(int argc, char **argv){

	float *mat, *vet, *vetres;
	int i, j;

	int id, nproc, resto = 0, ult_linha;

	MPI_Init(-np, 2);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);

	/* Aloca a matriz e os vetores */

	vet = (float *) malloc (N*sizeof(float));

	if ( id == 0 ){
		mat = (float *) malloc (N*N*sizeof(float));
		vetres = (float *)malloc(N*sizeof(float));
	}
	else{
		mat = (float *) malloc (N/nproc*N*sizeof(float));
		vetres = (float *) malloc (N/nproc*sizeof(float));
	}
	

	if ( id == 0 ){
	
		/* Inicializa a matriz */
		srand(time(NULL));	

		for( i=0 ; i<N ; i++ ){
			for ( j=0 ; j<N ; j++ ){
				mat[i*N+j] = i; // rand()%100;
			}
		}

		/* Inicializa o vetor */

		for( i=0 ; i<N ; i++ ){
			vet[i] = 1; //rand()%100;
		}

	}

	/* Distribui o vetor  */

	MPI_Bcast(vet,N,MPI_FLOAT,0,MPI_COMM_WORLD);

	/* Particiona a matriz  */

        if (id==0){
		resto = N % nproc;	
	}

	MPI_Scatter(&mat[resto*N], N/nproc*N, MPI_FLOAT, &mat[resto*N], N/nproc*N, MPI_FLOAT, 0, MPI_COMM_WORLD);

	/* Multiplicacao matriz x vetor  */

	ult_linha = N/nproc;
	ult_linha += resto;


	for( i=0 ; i < ult_linha ; i++ ){
		vetres[i] = 0;
		for ( j=0 ; j<N ; j++){
			vetres[i] += mat[i*N+j]*vet[j];
		}
	}
	
	/* Agrupar o vetor de resultados */

	MPI_Gather(&vetres[resto], N/nproc, MPI_FLOAT, &vetres[resto], N/nproc, MPI_FLOAT, 0, MPI_COMM_WORLD);


	/* Escreve o vetor */

	if ( id == 0 ){
		for( i=0 ; i<N ; i++ ){
			printf("vet[%d] = %f \n", i, vetres[i]);
		}
	}

	MPI_Finalize();

	
}

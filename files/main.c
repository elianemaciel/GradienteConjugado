#include "../include/entrada.h"
#include "../include/generate.h"
#include "../include/gradiente.h"

int main (int argc, char *argv[]) {

    unsigned int N;
    double *b, *Matriz, *x, tol, contIt;
    int nBandas, maxIter;
    char *caminho;
    FILE *arqSaida;

    if (Entrada(argc,argv,&N,&nBandas,&maxIter,&tol,&caminho) ) {

		double matrizSaida[maxIter][2];

		arqSaida = fopen(caminho,"w+");

		b = (double *) malloc (sizeof(double)*N);
		Matriz = (double *) malloc (sizeof(double)*(N*(nBandas/2+1)));
		x = (double *) malloc (sizeof(double)*N);

		srand(20162);
    
		generateMatriz (N, nBandas, Matriz);
		generateVector (N, b);

		fprintf(arqSaida,"###########\n");

		contIt = GradienteConjugado(Matriz,b,N,nBandas,maxIter,tol,x,arqSaida,matrizSaida);

		fprintf(arqSaida,"#\n");
		fprintf(arqSaida,"# Norma Euclidiana do Resíduo e Erro aproximado\n");
		imprimeSaida(matrizSaida,contIt,arqSaida);

		fprintf(arqSaida,"###########\n");
		fprintf(arqSaida,"%d\n",N);

		imprimeVetor(N,x,arqSaida);

		fclose(arqSaida);
		free (b);
		free (x);
		free (Matriz);
		return 0;
	}
	
	else
		return -1;

}

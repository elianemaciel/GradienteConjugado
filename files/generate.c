#include "../include/generate.h"

//Gera o vetor de termos indepedentes b
inline void generateVector ( unsigned int N, double *b) {
    double x; 
	int i;

    for (i = 0; i < N; i++) {
        x = i*M_PI/N;
        b[i] = 4*M_PI*M_PI *(  sin(2*M_PI*x) + sin (2*M_PI*(M_PI-x) ) );
    }
}

//Função que gera diagonal
inline int generateRandomDiagonal( unsigned int N, unsigned int k, unsigned int nBandas, double *diag ) {
    if ( !diag || N < 3 || nBandas > N/2 || k < 0 || k > nBandas )
        return (-1);

    /* garante valor dominante para diagonal principal */
    double fator = (k == 0) ? ((double)(nBandas-1)) : (0.0);

    double invRandMax = 1.0 / (double)RAND_MAX;
int i=0;
    for (i=0; i < N-k; ++i) {
        diag[i] = fator + (double)rand() * invRandMax;
    }

    return (0);
}

//Função que gera a matriz
inline void generateMatriz( unsigned int N, unsigned int nBandas, double *matriz ) {
    double vetorAux[N];
    int AuxTam = N, indMatriz = 0, k = nBandas/2;
int i = 0;
int j, indVetorAux = 0;
    for (i = 0; i <= k; i++) {
        generateRandomDiagonal (N, i, nBandas, vetorAux);
        for (j = indMatriz, indVetorAux = 0;  indVetorAux < AuxTam; indMatriz+= k+1, indVetorAux++, j++) {
            matriz[indMatriz] = vetorAux[indVetorAux];
        }

        indMatriz=i+1;
        AuxTam--;
    }
}

//Função que dado o indice da matriz, retorna o indice do elemento no vetor
inline int indexa (int i, int j, int nBandas) {

    if (j >= i)
	   return (nBandas+1)*i + (j - i);

	return (nBandas+1)*j + (i - j);
}

//Imprime Vetor
inline void imprimeVetor (unsigned int N, double *b, FILE *arqSaida) {
    int i = 0;
	for (i = 0; i < N; i++) {
        fprintf(arqSaida,"%.14g ", b[i]);
    }
    
	fprintf(arqSaida,"\n");
}

//Imprime saída no arquivo
inline void imprimeSaida(double matrizSaida[][2], int contIt, FILE *arqSaida) {
    int i;
    
	for (i = 0; i < contIt; i++) {
        fprintf(arqSaida, "# i=%d: %.14g %.14g\n",i,matrizSaida[i][iNorma],matrizSaida[i][iErro]);
    }
}

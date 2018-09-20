#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define iNorma 0
#define iErro 1

void generateVector (unsigned int N, double *b);
int generateRandomDiagonal(unsigned int N, unsigned int k, unsigned int nBandas, double *diag);
void generateMatriz(unsigned int N, unsigned int nBandas, double *matriz);
int indexa(int i, int j, int nBandas);
void imprimeVetor(unsigned int N, double *b, FILE *arqSaida);
void imprimeSaida(double matrizSaida[][2], int contIt, FILE *arqSaida);

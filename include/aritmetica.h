#include <stdlib.h>
#include <stdio.h>
#include "generate.h"

void multiplicaMatriz_Vetor(double *matriz, double *v, double *z, unsigned int N, int nBandas);
double multiplicaVetor_Vetor(double *a, double *b, unsigned int N);
void multiplicaInteiro_Vetor(double mult, double *v, double *vetorAux, unsigned int N);
void somaVetor(double *a, double *b, double *vetorFinal, unsigned int N);
void subtraiVetor(double *a, double *b, double *vetorFinal, unsigned int N);
void copiaVetor(double *a, double *b, unsigned int N);

#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "aritmetica.h"

// Função do tempo e do Gradiente
double timestamp(void);
int GradienteConjugado (double *matriz, double *b, unsigned int N, int nBandas, int maxIter, double tol, double *x, FILE *arqSaida, double matrizSaida[][2]);

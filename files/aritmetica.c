#include "../include/aritmetica.h"

//Multiplica Matriz com vetor
inline void multiplicaMatriz_Vetor(double *matriz, double *v, double *z, unsigned int N, int nBandas) {
    int i,j, k = nBandas/2;
    double sum;
    
	for (i = 0; i < N; i++) {
        sum = 0;
        
		for (j = 0; j < N; j++) {
            if (abs(i-j) <= k)
                sum += (matriz[indexa(i,j,k)] * v[j]);
        }
        
		z[i] = sum;
    }
}

//Multiplica vetor com vetor e retorna o resultado
inline double multiplicaVetor_Vetor (double *a, double *b, unsigned int N) {
    double result = 0;
    int i;

    for (i = 0; i < N; i++) {
        result += a[i]*b[i];
    }
    
	return result;

}

//Muiltiplica um número pelo vetor
inline void multiplicaInteiro_Vetor (double mult, double *v, double *vetorAux, unsigned int N) {
	int i;

	for (i = 0; i < N; i++) {
		vetorAux[i] = mult*v[i];
	}
}

//Soma 2 vetores
inline void somaVetor (double *a, double *b, double *vetorFinal, unsigned int N) {
	int i;

	for (i = 0; i < N; i++) {
		vetorFinal[i] = a[i]+b[i];
	}
}

//Subtrai 2 vetores
inline void subtraiVetor (double *a, double *b, double *vetorFinal, unsigned int N) {
	int i;

	for (i = 0; i < N; i++) {
		vetorFinal[i] = a[i]-b[i];
	}
}

//Copia um vetor
inline void copiaVetor (double *a, double *b, unsigned int N) {
	int i;
	
	for (i = 0; i < N; i++) {
		b[i] = a[i];
	}

}


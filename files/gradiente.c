#include "../include/gradiente.h"

inline double timestamp(void) {
    struct timeval tp;
    gettimeofday(&tp, NULL);
    return((double)(tp.tv_sec*1000.0 + tp.tv_usec/1000.0));
}

//Metodo do Gradiente Conjugado
//Matriz = Vetor que simula a Matriz gerada
//B = Vetor de termos indepentes
//N = tamanho do sistema
//nBandas = Número de bandas presentes na matriz
//maxIter = Número maximo de iterações, caso não seja definida é igual a N
//tol = Número de tolerancia de erro permitida, caso não seja definida o erro e 0.00001
//x = Vetor solução do sistema
//arqSaida = Arquivo de Saida que sera guardado o tempo e o vetor solução
//matrizSaida = Matriz que guardara o Erro e a Norma de cada iteração
inline int GradienteConjugado (double *matriz, double *b, unsigned int N, int nBandas, int maxIter, double tol, double *x, FILE *arqSaida, double matrizSaida[][2]) {

    //aux = variavel auxiliar
	//erro = erro aproximado do sistema
	//s. v, z, vetorAux = Vetor Auxiliares
	//r = Vetor do residuo
	double aux,erro, m, s;
    double *r, *v, *z, *vetorAux;
    int k;

	//variaveis para medição de tempo
    double maxtmetodo, mintmetodo, conttmetodo, mediatmetodo;
    double tminicio, tmfinal;
    double maxtresd, mintresd, conttresd, mediatresd;
    double tresdinicio, tresdfinal;

    conttmetodo = 0.0, conttresd = 0.0;
    maxtmetodo = 0.0, maxtresd = 0.0;
    mintmetodo = 1000000, mintresd = 1000000;

    r = (double*) malloc (sizeof (double) * N);
    z = (double*) malloc (sizeof (double) * N);
    v = (double*) malloc (sizeof (double) * N);
    vetorAux = (double*) malloc (sizeof (double) * N);

    copiaVetor (b, r, N);
    copiaVetor (b, v, N);
    memset(x,0,sizeof(x)*N);


    aux = multiplicaVetor_Vetor (r,r, N);

    for (k = 0; k < maxIter; k++) {

        tminicio = timestamp();

        multiplicaMatriz_Vetor (matriz, v, z, N, nBandas);

        s = aux/(multiplicaVetor_Vetor(v,z, N));

        multiplicaInteiro_Vetor (s, v, vetorAux, N);
        somaVetor(vetorAux,x,x,N);

        multiplicaInteiro_Vetor (s, z, vetorAux, N);
        subtraiVetor (r, vetorAux, r, N);

        tresdinicio = timestamp();
        
		erro = multiplicaVetor_Vetor (r, r, N);

        tresdfinal = timestamp() - tresdinicio;
        tmfinal = (timestamp() - tminicio);

        conttmetodo += tmfinal;
        conttresd += tresdfinal;

        if ( tmfinal < mintmetodo)
            mintmetodo = tmfinal;

        if ( tmfinal > maxtmetodo)
            maxtmetodo = tmfinal;

        if ( tresdfinal > maxtresd )
        	maxtresd = tresdfinal;

        if ( tresdfinal < mintresd)
            mintresd = tresdfinal;

        matrizSaida[k][iNorma] = sqrt(erro);
        matrizSaida[k][iErro] = erro;
    
        if ( erro <= tol) {
            k = k+1;
            break;
        }

        m = erro/aux;
        aux = erro;

        multiplicaInteiro_Vetor(m,v, vetorAux, N);
        somaVetor (r, vetorAux, v, N);

    }

    mediatmetodo = conttmetodo/k;
    mediatresd = conttresd/k;

    fprintf(arqSaida, "# Tempo Método CG: %lf %lf %lf\n", mintmetodo, mediatmetodo, maxtmetodo);
    fprintf(arqSaida, "# Tempo Resíduo: %lf %lf %lf\n", mintresd, mediatresd, maxtresd);

    free (z);
    free (r);
    free (v);
    free (vetorAux);
    return k;
}

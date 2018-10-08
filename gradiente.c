#include <sys/time.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "hb_io.h"
#include <math.h>

#define iNorma 0
#define iErro 1

//Função que dado o indice da matriz, retorna o indice do elemento no vetor
inline int indexa (int i, int j, int nBandas) {

    if (j >= i)
	   return (nBandas+1)*i + (j - i);

	return (nBandas+1)*j + (i - j);
}

//Multiplica Matriz com vetor
inline void multiplicaMatriz_Vetor(int *matriz, double *v, double *z, unsigned int N, int nBandas) {
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
inline void copiaVetor (int *a, double *b, unsigned int N) {
	int i;

	for (i = 0; i < N; i++) {
		b[i] = a[i];
	}

}


inline int Entrada(int argc, char *argv[], unsigned int *N, int *nBandas, int *maxIter, double *tol) {
    //lendo argumento n
    if (argv[1] != NULL)
        *N = atoi (argv[1]);
    else {
        fprintf (stderr, "ERRO, argumento N e nBandas não informados\n");
		return 0;
	}

    //lendo argumento k
    if (argv[2] != NULL)
        *nBandas = atoi (argv[2]);
    else {
        fprintf (stderr, "ERRO, argumento nBandas não informado\n");
		return 0;
	}


    //caso -i exista, le argumento maxiter, se nao faz maxiter = n
    if ( (argv[3] != NULL) && !strcmp(argv[3], "-i") )
        *maxIter = atoi (argv[4]);
    else
        *maxIter = *N;

    //caso -i n exista, ve se terceiro argumento eh -t
    if ( (argv[3] != NULL) && !strcmp(argv[3], "-t") )
        *tol = atof (argv[4]);

    //caso -i exista, ve se o quinto argumento eh -t
    else if ( (argv[5] != NULL) && !strcmp (argv[5], "-t") )
        *tol = atof (argv[6]);

    //tolerancia n foi definid
    else
        *tol = 4.9e-324;

	return 1;
}

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

    int i = 0;
    for (i=0; i < N-k; ++i) {
        diag[i] = fator + (double)rand() * invRandMax;
    }

    return (0);
}

//Função que gera a matriz
inline void generate_matriz( unsigned int N, unsigned int nBandas, double *matriz ) {
    double vetorAux[N];
    int AuxTam = N, indMatriz = 0, k = nBandas/2;
    int i = 0, j, indVetorAux = 0;
    for (i = 0; i <= k; i++) {
        generateRandomDiagonal (N, i, nBandas, vetorAux);

        for (j = indMatriz, indVetorAux = 0;  indVetorAux < AuxTam; indMatriz+= k+1, indVetorAux++, j++) {
            matriz[indMatriz] = vetorAux[indVetorAux];
        }

        indMatriz=i+1;
        AuxTam--;
    }
}


//Imprime Vetor
inline void imprimeVetor (unsigned int N, double *b) {
      int i = 0;
	for (i = 0; i < N; i++) {
        printf("%.14g ", b[i]);
    }

	printf("\n");
}

//Imprime saída no arquivo
inline void imprimeSaida(double matrizSaida[][2], int contIt) {
    int i;

	for (i = 0; i < contIt; i++) {
        printf("# i=%d: %.14g %.14g\n",i,matrizSaida[i][iNorma],matrizSaida[i][iErro]);
    }
}


//Metodo do Gradiente Conjugado
//Matriz = Vetor que simula a Matriz gerada
//B = Vetor de termos indepentes
// N = tamanho do sistema
//nBandas = Número de bandas presentes na matriz
//maxIter = Número maximo de iterações, caso não seja definida é igual a N
//tol = Número de tolerancia de erro permitida, caso não seja definida o erro e 0.00001
//x = Vetor solução do sistema
//matrizSaida = Matriz que guardara o Erro e a Norma de cada iteração
inline int GradienteConjugado (int *matriz, int *b, unsigned int N, int nBandas, int maxIter, double tol, double *x, double matrizSaida[][2]) {

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

        // tminicio = timestamp();

        multiplicaMatriz_Vetor (matriz, v, z, N, nBandas);

        s = aux/(multiplicaVetor_Vetor(v,z, N));

        multiplicaInteiro_Vetor (s, v, vetorAux, N);
        somaVetor(vetorAux,x,x,N);

        multiplicaInteiro_Vetor (s, z, vetorAux, N);
        subtraiVetor (r, vetorAux, r, N);

        // tresdinicio = timestamp();

		erro = multiplicaVetor_Vetor (r, r, N);

        // tresdfinal = timestamp() - tresdinicio;
        // tmfinal = (timestamp() - tminicio);

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

    printf("# Tempo Método CG: %lf %lf %lf\n", mintmetodo, mediatmetodo, maxtmetodo);
    printf("# Tempo Resíduo: %lf %lf %lf\n", mintresd, mediatresd, maxtresd);

    free (z);
    free (r);
    free (v);
    free (vetorAux);
    return k;
}

int main (int argc, char *argv[]) {
    unsigned int N;
    double *b, *Matriz, *x, tol;
    int nBandas, maxIter, contIt;

    maxIter = 100;
    tol = 0.2;
    // Nova lib
    int *colptr = NULL;
    int indcrd;
    char *indfmt = NULL;
    FILE *input;
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


    input = fopen("matrizMenor.rsa", "rt" );

    if ( !input ){
      printf ( "\n" );
      printf ( "  Error opening the file.\n" );
    }

    hb_header_read(input, &title, &key, &totcrd, &ptrcrd, &indcrd,
        &valcrd, &rhscrd, &mxtype, &nrow, &ncol, &nnzero, &neltvl, &ptrfmt,
        &indfmt, &valfmt, &rhsfmt, &rhstyp, &nrhs, &nrhsix
    );

    colptr = ( int * ) malloc ( ( ncol + 1 ) * sizeof ( int ) );

    if ( mxtype[2] == 'A' )
    {
      rowind = ( int * ) malloc ( nnzero * sizeof ( int ) );
    }
    else if ( mxtype[2] == 'E' )
    {
      rowind = ( int * ) malloc ( neltvl * sizeof ( int ) );
    }
    else
    {
      printf ( "\n" );
      printf ( "TEST05 - Warning!\n" );
      printf ( "  Illegal value of MXTYPE character 3.\n" );
    }


    hb_structure_read ( input, ncol, mxtype, nnzero, neltvl,
      ptrcrd, ptrfmt, indcrd, indfmt, colptr, rowind );

    if ( mxtype[2] == 'A' )
    {
      values = ( double * ) malloc ( nnzero * sizeof ( double ) );
    }
    else if ( mxtype[2] == 'E' )
    {
      values =  ( double * ) malloc ( neltvl * sizeof ( double ) );
    }
    else
    {
      printf ( "\n" );
      printf ( "  Erro = '%c'\n", mxtype[2] );
    }

    hb_values_read ( input, valcrd, mxtype, nnzero, neltvl, valfmt, values );

    fclose ( input );

    printf ( "\n" );
    printf ( "  '%s'\n", title );
    printf ( "  KEY =    '%s'\n", key );
    printf ( "\n" );
    printf ( "  NROW =   %d\n", nrow );
    printf ( "  NCOL =   %d\n", ncol );
    printf ( "  NNZERO = %d\n", nnzero );
    printf ( "  NELTVL = %d\n", neltvl );

    int i=0;

    // Vetor colptr
	for(i=0;i<=nrow;i++){
		printf ( "  colptr =   %d\n", colptr[i] );
	}
	i=0;
	printf("\n");

    // Matriz rowind
	for(i=0; i<nnzero;i++){
		printf ( "  rowind =   %d\n", rowind[i] );
	}

	i=0;

    // Valores
	printf("\n");
    for(i=0; i<nnzero;i++){
		printf ( "  values =   %.2f\n", values[i] );
	}

    /*
     Print out the  header information.
    */
    // hb_header_print ( title, key, totcrd, ptrcrd, indcrd, valcrd,
    //   rhscrd, mxtype, nrow, ncol, nnzero, neltvl, ptrfmt, indfmt, valfmt,
    //   rhsfmt, rhstyp, nrhs, nrhsix );


	double matrizSaida[maxIter][2];


	// contIt = GradienteConjugado(rowind,colptr,nnzero,nnzero,maxIter,tol,values,matrizSaida);

	// printf("#\n");
	// printf("# Norma Euclidiana do Resíduo e Erro aproximado\n");
	//imprimeSaida(matrizSaida,contIt);

    free ( colptr );
    free ( rowind );
    free ( values );

	return 0;

}

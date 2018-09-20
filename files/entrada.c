#include "../include/entrada.h"

inline int Entrada(int argc, char *argv[], unsigned int *N, int *nBandas, int *maxIter, double *tol, char **caminho) {
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


    //caso -i e -t n existam
    if ( (argv[3] != NULL) && !strcmp (argv[3], "-o") )
        *caminho = argv[4];

    //caso -i ou -t n existam
    else if ( (argv[5] != NULL) && !strcmp (argv[5], "-o") )
        *caminho = argv[6];

    //caso -i e -t existam
    else if ( (argv[7] != NULL) && !strcmp (argv[7], "-o") )
        *caminho = argv[8];
    else {
        fprintf (stderr, "ERRO, caminho de saida não informado\n");
		return 0;
	}
	return 1;
}

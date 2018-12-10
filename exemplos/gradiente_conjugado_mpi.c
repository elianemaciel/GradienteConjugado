/*###########################################################################
	AUTOR:	MARCAL NUNES DE OLIVEIRA
			RENAN FARAON
############################################################################*/

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define MM_MAX_LINE_LENGTH 1025
#define MM_MAX_TOKEN_LENGTH 64

/********************* Matrix Market error codes ***************************/


#define MM_COULD_NOT_READ_FILE	11
#define MM_PREMATURE_EOF		12
#define MM_NOT_MTX				13
#define MM_NO_HEADER			14
#define MM_UNSUPPORTED_TYPE		15
#define MM_LINE_TOO_LONG		16
#define MM_COULD_NOT_WRITE_FILE	17

/***************************************************************************/

typedef char MM_typecode[4];

//le uma linha e não faz nada 
int mm_read_banner(FILE *f)
{
	char line[MM_MAX_LINE_LENGTH];
	if (fgets(line, MM_MAX_LINE_LENGTH, f) == NULL) 
        return MM_PREMATURE_EOF;
}

/****************************************************************************
				le os parametros da Matriz
C = Total de Linhas da Matriz
E = Número de elementos não nulos
*****************************************************************************/
int mm_read_mtx_size(FILE *f, int *C, int *E, int *V)
{
	char line[MM_MAX_LINE_LENGTH];
	char s[3];
	if (fgets(line, MM_MAX_LINE_LENGTH, f) == NULL) 
        return MM_PREMATURE_EOF;
    
    sscanf(line, "%s %d %*s %d %d",s,  C, E, V);
    
    if(s[1] != 'S'){
    	printf("Formato de arquivo não reconhecido");
    	exit(1);    	
	}
    
    return 0;
}

/****************************************************************************
				le os primeiros parametros do arquivo
M = Total de Linhas do arquivo (sem o cabeçalho)
N = Número de linhas do vetor row_ind
nz = Número de linhas do vetor col_ptr
K = Número de linhas do vetor value
*****************************************************************************/
int mm_read_size_row(FILE *f, int *M, int *N, int *nz , int *k)
{
    char line[MM_MAX_LINE_LENGTH];
    int num_items_read;

    /* set return null parameter values, in case we exit with errors */
    *M = *N = *nz = *k = 0;

    if (fgets(line,MM_MAX_LINE_LENGTH,f) == NULL) 
        return MM_PREMATURE_EOF;

    /* line[] is either blank or has M,N, nz, k */
    sscanf(line, "%d %d %d %d", M, N, nz,k);        

    return 0;
}

/******************************************************************
			le as posicoes das linhas
******************************************************************/			
int mm_read_row(FILE *f, int *row_ind, int E){
	char line[MM_MAX_LINE_LENGTH];
	int i;
	for(i = 0; i <E; i++){		      	
		fscanf(f, "%d", &row_ind[i]);
		row_ind[i]--;
	}	
	return 0;
}

/******************************************************************
			le as posicoes das colunas
******************************************************************/
int mm_read_column(FILE *f, int *col_ptr, int n){
	char line[MM_MAX_LINE_LENGTH];
	int i,k;
	for(i = 0; i <=n; i++){	
  		fscanf(f, "%d", &col_ptr[i]);
  		col_ptr[i] --;
	}	
	return 0;
}

/******************************************************************
			le os valores da matriz
******************************************************************/
int mm_read_mtx(FILE *f, double *val, int E){
	char line[MM_MAX_LINE_LENGTH];
	int i;
	for(i = 0; i <E; i++){		      	
		fscanf(f, "%lf", &val[i]);
	}	
	return 0;
}

/******************************************************************
			le os valores do vetor
******************************************************************/
int mm_read_vet(FILE *f, double *vet, int n){
	char line[MM_MAX_LINE_LENGTH];
	int i;
	for(i = 0; i <n; i++){		      	
		fscanf(f, "%lf", &vet[i]);
	}	
	return 0;
}

/*-------------------------------------------------------------------------------------------------------*/
float multiplicacao_vetor_vetor(double *vetor1, double *vetor2, int n){
	int i;
	float result = 0;
	for (i = 0; i < n; i++) 
	{ 
	    result+= vetor1[i] * vetor2[i];
	}
	return result;
}


/*-------------------------------------------------------------------------------------------------------*/
void multiplicacao_matriz_vetor (double *val, int *col_ptr, int *row_ind, double *vet, double *res, int n, int nc, int id){
	int i,j,l;
	//verificar a possibilidade de gerar um ponteiro novo
	for(i=0;i<id;i++)
		res[i] = 0;
		
	for(i=n,l=0;i<nc;i++,l++){
		for(j=col_ptr[l];j<col_ptr[l+1];j++){
			res[i] += val[j] * vet[row_ind[j]];
			if( i  != row_ind[j]){
				res[row_ind[j]] += val[j] * vet[i];
			}
		}
	}		
}



/*-------------------------------------------------------------------------------------------------------*/
void copiaVetor(double *vetor, double *copia, int n){
	int i;
	for (i = 0; i < n; i++) 
	{ 
	    copia[i] = vetor[i];
	}
}

/*-------------------------------------------------------------------------------------------------------*/
void multiplicacao_vetor_double(double *vetor1, double num, double *res, int n){
	int i;
	for(i = 0 ; i < n ; i++){
		res[i] = vetor1[i] * num;
	}
	
}

/*-------------------------------------------------------------------------------------------------------*/
void soma_vetor_double(double *vetor1, double num, double *res, int n){
	int i;
	for(i = 0 ; i < n ; i++){
		res[i] = vetor1[i] + num;
	}
	
}

/*-------------------------------------------------------------------------------------------------------*/
void soma_vetor_vetor(double *vetor1, double *vetor2, double *res, int n){
	int i;
	for(i = 0 ; i < n ; i++){
		res[i] = vetor1[i] + vetor2[i];
	}
	
}

/*-------------------------------------------------------------------------------------------------------*/
void subtracao_vetor_vetor(double *vetor1, double *vetor2, double *res, int n){
	int i;
	for(i = 0 ; i < n ; i++){
		res[i] = vetor1[i] - vetor2[i];
	}
	
}


int main(int argc, char *argv[]){

	int 	ret_code;
	MM_typecode matcode;
	FILE 	*f;
	int 	M, N, nz, K, E, V;      
	int 	j,l;
    int 	imax 		= 100000;
	float 	erro 		= 0.00001; 
	int id, np;
	MPI_Status s;
	
	
	int 	i;
	int *row_ind;
	int *col_ptr;
	double 	*val;

	int 	n;
	double 	*r 			= NULL;
	double 	*res		= NULL;
	double 	*res2		= NULL;
	double 	*b 			= NULL;
	double 	sigma_novo;
	double 	sigma_velho;
	double 	sigma0;
	double 	*d 			= NULL;
	double 	*q 			= NULL;
	double 	alpha ;
	double  *x;
	


	if (argc < 1){
		fprintf(stderr, "%s < arquivo matriz >\n", argv[0]);
		exit(1);
	}
    
	f = fopen(argv[1], "r");

	if ( f == NULL ){
		printf("Erro ao abrir o arquivo\n");
		exit(1);
	}
	
    /**************************************************
    					MPI
    **************************************************/
    MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &np);

	if(id == 0){
		// Le os Banner do arquivo
		mm_read_banner(f);
	
		// Le a primeira linha dos parametros da matriz
		mm_read_size_row(f, &M, &N, &nz, &K);
		// Le a segunda linha dos parametros da matriz
		mm_read_mtx_size(f,&n,&E,&V);
		
		
	}
	
	//envia para todos processos os tamanhos da matriz e vetor
	MPI_Bcast(&E, 1, MPI_INT, 0, MPI_COMM_WORLD);	
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	//divide o tamanho do vetor pelos processos
	int div= (1+n)/np;
		

	if(id == 0){//principal
		// Aloca memoria para armazenar a matriz
		row_ind = (int *) malloc(sizeof(int) * E);
		val = (double *) malloc(sizeof(double) * E); 		
	    col_ptr = (int *) malloc(sizeof(int) * (n+1));
	    b = (double *) malloc(sizeof(double) * n); 
	}
	else{
		//somente o col_ptr, pois os demais serão variaveis conforme o col_ptr receber os valores
		col_ptr = (int *) malloc(sizeof(int) * (div+1));
	}
		
	//b = (double *) malloc(sizeof(double) * n);  
	res = (double *) malloc(sizeof(double) * n);
	
	
	if(id ==0){	
		//pula uma linha
		mm_read_banner(f);
		
		//verifica se existe o vetor no arquivo, caso contrario gera com 1
		if(V){
			// Le os valores do vetor
	    	mm_read_vet(f, b, n );
		}
		else{
			for (i = 0; i < n; i++) {
		        b[i] = 1;    
		    }
		}
		printf("lendo colunas\n");
	    // le os indicadores de onde começa e termina cada linha dentro do vetor de linhas
	    mm_read_column(f, col_ptr, n );
	    	
		printf("lendo linhas %d\n",E);
	    // Le as linhas
	    mm_read_row(f, row_ind, E );
	    printf("lendo valores\n");
	    //le os valores da matriz
	    mm_read_mtx(f, val, E );
		
	 	fclose(f);
	 	
	 	r = (double *) malloc(sizeof(double ) * n);
		q = (double *) malloc(sizeof(double ) * n);
		res2 = (double *) malloc(sizeof(double) * n);

	}
	
	int *displs, *scounts, *sdiv, *sdivC;
	displs = (int *)malloc(np*sizeof(int)); 
	scounts = (int *)malloc(np*sizeof(int)); 
	/**************************************************************************
		cria os vetores com as informações
		displs - vetor com o valor de onde começa cada processo dentro do col_ptr
		scounts - vetor com a quantidade que será enviada para cada processo
	**************************************************************************/
	for (i=0; i<np; ++i) { 
		displs[i] = i*(div); 
		scounts[i] = div+1; 
	} 

	if(id == 0)
		printf("enviando col_ptr\n");	
	//envia o vetor dividido utilizando um offset
	MPI_Scatterv(col_ptr,scounts, displs,MPI_INT,col_ptr,div+1,MPI_INT,0,MPI_COMM_WORLD);

	/***********************************************************************************
	Para o servidor:
		sdiv - tem o valor do col_ptr que inicia  o processo i 
		sdivC - tem a quantidade que vai ser enviada para os clientes, porem o cálculo 
				da ultima é feito de forma diferente pois pode ter menos o col_ptr

	Para Cliente:
		sdiv - tem o valor do col_ptr que inicia o processo, sempre busca do 0, pois como foi
				dividido cada processo vai sempre iniciar no col_ptr[0]
		sdivC - tem a quantidade que vai ser recebida, porem o cálulo da última também vai
				ser diferente pois vai ter o delimitador o E.
	***********************************************************************************/
	if(id == 0 ){		
		sdiv = (int *)malloc(np*sizeof(int)); 
		sdivC = (int *)malloc(np*sizeof(int)); 
		for (i=0; i<np; ++i) { 
			sdiv[i] = col_ptr[i*(div)];
			if(i != np-1)
				sdivC[i] = col_ptr[(i*(div))+(div)] - sdiv[i];
			else
				sdivC[i] = col_ptr[n] - sdiv[i];
		} 
	}
	else{		
		sdiv = (int *)malloc(1*sizeof(int)); 
		sdivC = (int *)malloc(1*sizeof(int));
		sdiv[0] = col_ptr[0];
		if(id == np-1)
			sdivC[0] =  E - col_ptr[0];
		else
			sdivC[0] = col_ptr[div] - col_ptr[0];
			
		row_ind = (int *) malloc(sizeof(int) * sdivC[0]);
		val = (double *) malloc(sizeof(double) * sdivC[0]); 	
		int ini = col_ptr[0];
		//alterado o valor do col_ptr, pois no row_ind sempre vai começar no 0 para cada processo
		for(i=0; i<=div;i++)
			col_ptr[i] = col_ptr[i] - ini;
	}
	
	//Envia o row_ind
	if(id == 0){
		printf("enviando row_ind\n");
		for(i=1;i<np;i++)
			MPI_Send(&row_ind[sdiv[i]],sdivC[i],MPI_INT, i,100,MPI_COMM_WORLD);
	}
	else{
		MPI_Recv(row_ind,sdivC[0],MPI_INT, 0,100,MPI_COMM_WORLD,&s);
	}

	//envia o indicador de limite de linha
	int row_limit = 0;
	if(id == 0){
		row_limit = displs[1];
		for(i=1;i<np-1;i++)
			MPI_Send(&displs[i+1],1,MPI_INT, i,300,MPI_COMM_WORLD);
		displs[np-1]++;
		MPI_Send(&displs[np-1],1,MPI_INT, i,300,MPI_COMM_WORLD);
	}
	else{
		MPI_Recv(&row_limit,1,MPI_INT, 0,300,MPI_COMM_WORLD,&s);
	}
	

	//Envia o vetor val
	if(id == 0){
		printf("enviando valores\n");
		for(i=1;i<np;i++)
			MPI_Send(&val[sdiv[i]],sdivC[i],MPI_DOUBLE, i,200,MPI_COMM_WORLD);
	}
	else{
		MPI_Recv(val,sdivC[0],MPI_DOUBLE, 0,200,MPI_COMM_WORLD,&s);
	}

	
	// x = zeros(n,1);
	x = (double *) malloc(sizeof(double ) * n);
	for (i = 0; i < n; i++) {
		x[i] = 0;    
    	}
	d = (double *) malloc(sizeof(double ) * n);
	

	//r = b - A * x;
	//res = A * x
	multiplicacao_matriz_vetor(val, col_ptr, row_ind, x, res, displs[id], row_limit,n);

	if(id == 0){
		for(i=1;i<np;i++){
			MPI_Recv(res2,n,MPI_DOUBLE, i,10,MPI_COMM_WORLD,&s);
			soma_vetor_vetor(res,res2,res,n);
		}
		//r = b - res
		subtracao_vetor_vetor(b, res, r, n);

		//d=r
		copiaVetor(r,d,n);
		
		//sigma_novo = r' * r;
		sigma_novo = multiplicacao_vetor_vetor(r, r,n);
		sigma0 = sigma_novo;
		
		
		l = 1;
		j = 1;	
		while(j < imax && sigma_novo > ((erro*erro)*sigma0)){
			for(i=1;i<np;i++){
				MPI_Send(&l,1,MPI_INT, i,20,MPI_COMM_WORLD);
				MPI_Send(d,n,MPI_DOUBLE, i,200,MPI_COMM_WORLD);	
			}

			// q = A * d;
			multiplicacao_matriz_vetor(val, col_ptr, row_ind, d, q, displs[id], row_limit,n);
			for(i=1;i<np;i++){
				MPI_Recv(res2,n,MPI_DOUBLE, i,10,MPI_COMM_WORLD,&s);
				soma_vetor_vetor(q,res2,q,n);
			}
			//alpha = sigma_novo/(d' * q); 
			alpha = sigma_novo/multiplicacao_vetor_vetor(d, q,n);
			
			//x = x + alpha * d;
			// res = alpha * d			
			multiplicacao_vetor_double( d, alpha, res, n);
			// x = x + res;
			soma_vetor_vetor(x, res, x, n);
			
			if(j%50 == 0){
				for(i=1;i<np;i++){
					MPI_Send(&l,1,MPI_INT, i,20,MPI_COMM_WORLD);				
					MPI_Send(x,n,MPI_DOUBLE, i,200,MPI_COMM_WORLD);				
				}
				//r = b - A * x;
				//res = A * x
				multiplicacao_matriz_vetor(val, col_ptr, row_ind, x, res, displs[id], row_limit,n);
				for(i=1;i<np;i++){
					MPI_Recv(res2,n,MPI_DOUBLE, i,10,MPI_COMM_WORLD,&s);
					soma_vetor_vetor(res,res2,res,n);
				}
				//r = b - res
				subtracao_vetor_vetor(b, res, r, n);
			}
			else{
				//r = r - alpha * q;
				// res = alpha * q			
				multiplicacao_vetor_double( q, alpha, res, n);
				// r = r - res
				subtracao_vetor_vetor(r, res, r,  n);
			}
			 
			sigma_velho = sigma_novo;
			sigma_novo = multiplicacao_vetor_vetor(r, r,n);
			double beta = sigma_novo / sigma_velho;
			
			// d = r + beta * d;
			// res = beta * d
			multiplicacao_vetor_double(d, beta,res, n);
			// d = r + res;
			soma_vetor_vetor(r, res, d, n);

				
			j++;
		}
		// envio de 1 para informar aos processos que terminou.
		l=0;
		for(i=1;i<np;i++)
			MPI_Send(&l,1,MPI_INT,i,20,MPI_COMM_WORLD);

		 
		printf("Finalizou com o i = %d\n\n", j);
		printf("Vetor x =\n");
		for (i = 0; i < n; i++) 
		 { 
			printf("%lf\n", x[i]);    
		 }
	}
	else{
		l = 1;
		MPI_Send(res,n,MPI_DOUBLE, 0,10,MPI_COMM_WORLD);
		MPI_Recv(&l,1,MPI_INT, 0,20,MPI_COMM_WORLD,&s);
		while(l){
			MPI_Recv(x,n,MPI_DOUBLE, 0,200,MPI_COMM_WORLD,&s);
			multiplicacao_matriz_vetor(val, col_ptr, row_ind, x, res, displs[id], row_limit,n);
			MPI_Send(res,n,MPI_DOUBLE, 0,10,MPI_COMM_WORLD);
			MPI_Recv(&l,1,MPI_INT, 0,20,MPI_COMM_WORLD,&s);
		}
	}

	//desalocar
	free(row_ind);
	free(col_ptr);
	free(val);
	free(r);
	free(b);
	free(d);
	free(q);
	free(x);
	free(res);
	free(res2);

	MPI_Finalize();
}

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


#include "global_functions_MULTI_KP.h"

#define print_ist_features
//#define print_ist_features_CONFLICTS
//#define print_sol
//#define print_sol_frac

// return a integer random value in range min-max
/*****************************************************************/
int randomBETWEEN(int min,int max)
/*****************************************************************/
{
	return rand() % (max - min +1) +min;
}

// return random value in range min-max
/*****************************************************************/
double randomBETWEEN_double(int min,int max)
/*****************************************************************/
{
	return (rand()/(double)RAND_MAX)*(max-min) + min;
}



// return a random value in range 0.0-1.0
/*****************************************************************/
double random01()
/*****************************************************************/
{
	return ((double) rand() / RAND_MAX);
}

/*********************************/
int read_instance_file_multi_kp(data_MULTI_KP *MULTI_KP_instance,double conflict_density)
/*********************************/
{


	ifstream in(MULTI_KP_instance->instname);
	if(!in)
	{
		cout << "File could not be opened. " << endl;
		exit(1);
	}


	int dummy;

	in >> MULTI_KP_instance->N_ITEMS;
	in >> MULTI_KP_instance->N_KNAPSACKS;
	in >> dummy;



	cout << "\n\nFeatures:\n";
	cout << "N_ITEMS\t\t" << MULTI_KP_instance->N_ITEMS << "\n";
	cout << "N_KNAPSACKS\t" << MULTI_KP_instance->N_KNAPSACKS << "\n";

	MULTI_KP_instance->MULTI_PROFITS=(double*)malloc(sizeof(double)*MULTI_KP_instance->N_ITEMS);


	for(int i=0;i<MULTI_KP_instance->N_ITEMS;i++)
	{

		in >> MULTI_KP_instance->MULTI_PROFITS[i];

	}


	MULTI_KP_instance->MATRIX_WEIGHTS=(int**)malloc(sizeof(int*)*MULTI_KP_instance->N_KNAPSACKS);

	for(int k=0;k<MULTI_KP_instance->N_KNAPSACKS;k++)
	{
		MULTI_KP_instance->MATRIX_WEIGHTS[k]=(int*)malloc(sizeof(int)*MULTI_KP_instance->N_ITEMS);
	}

	for(int j=0;j<MULTI_KP_instance->N_KNAPSACKS;j++)
	{
		for(int i=0;i<MULTI_KP_instance->N_ITEMS;i++)
		{
			in >> MULTI_KP_instance->MATRIX_WEIGHTS[j][i];
		}
	}

	MULTI_KP_instance->RHS_CAPACITIES=(int*)malloc(sizeof(int)*MULTI_KP_instance->N_KNAPSACKS);

	for(int j=0;j<MULTI_KP_instance->N_KNAPSACKS;j++)
	{
		in >> MULTI_KP_instance->RHS_CAPACITIES[j];
	}


	////////////////////////////////////////////////////////////////////////////////////////////////////
	int n_conflicts=0;

	MULTI_KP_instance->CONF_MATRIX=new int*[MULTI_KP_instance->N_ITEMS];
	for ( int j = 0; j < MULTI_KP_instance->N_ITEMS; j++ )
	{
		MULTI_KP_instance->CONF_MATRIX[j]=new int[MULTI_KP_instance->N_ITEMS];
	}

	for ( int j = 0; j < MULTI_KP_instance->N_ITEMS; j++ )
	{
		for ( int i = 0; i < MULTI_KP_instance->N_ITEMS; i++ )
		{
			MULTI_KP_instance->CONF_MATRIX[j][i]=0;
		}
	}

	for(int i=0; i<MULTI_KP_instance->N_ITEMS; i++)
	{

		for(int j=i+1; j<MULTI_KP_instance->N_ITEMS; j++)
		{

			double n_random=randomBETWEEN_double(0,1);

			if(n_random<=conflict_density)
			{

				//					cout << "n_random\t" << n_random << endl;
				//					cout << "conflict_density\t" << conflict_density << endl;
				//					cout << "conflict\t" << i << "\t" << j << endl;

				MULTI_KP_instance->CONF_MATRIX[i][j]=1;
				MULTI_KP_instance->CONF_MATRIX[j][i]=1;

				n_conflicts++;
			}
		}
	}
	////////////////////////////////////////////////////////////////////////////////////////////////////



#ifdef print_ist_features

	cout << "MULTI_PROFITS\n";
	for(int i=0;i<MULTI_KP_instance->N_ITEMS;i++)
	{
		cout <<  MULTI_KP_instance->MULTI_PROFITS[i] << "\t";
	}
	cout << endl;

	cout << "RHS_CAPACITIES\n";
	for(int j=0;j<MULTI_KP_instance->N_KNAPSACKS;j++)
	{
		cout <<  MULTI_KP_instance->RHS_CAPACITIES[j] << "\t";
	}
	cout << endl;

	cout << "MATRIX_WEIGHTS\n";
	for(int j=0;j<MULTI_KP_instance->N_KNAPSACKS;j++)
	{
		for(int i=0;i<MULTI_KP_instance->N_ITEMS;i++)
		{
			cout <<  MULTI_KP_instance->MATRIX_WEIGHTS[j][i] << "\t";
		}
		cout << endl;
	}
	cout << endl;


#ifdef	print_ist_features_CONFLICTS
	cout << "CONFLICTS\n\n";
	for ( int j = 0; j < MULTI_KP_instance->N_ITEMS; j++ )
	{
		for ( int i = 0; i < MULTI_KP_instance->N_ITEMS; i++ )
		{
			cout << MULTI_KP_instance->CONF_MATRIX[j][i] << " ";
		}
		cout << "\n";
	}
	cout << "\n";
#endif

#endif

	in.close();



	cout << "Instance read!\n";
	///////////////////////////////////////////////////////////


	return n_conflicts;
}


/*********************************/
void generate_instance_MULTI_KP(data_MULTI_KP *MULTI_KP_instance)
/*********************************/
{

	MULTI_KP_instance->N_ITEMS=MULTI_KP_instance->number_of_items;
	MULTI_KP_instance->N_KNAPSACKS=MULTI_KP_instance->number_of_KPs;


	MULTI_KP_instance->MULTI_PROFITS=(double*)malloc(sizeof(double)*MULTI_KP_instance->N_ITEMS);

	MULTI_KP_instance->MATRIX_WEIGHTS=(int**)malloc(sizeof(int*)*MULTI_KP_instance->N_KNAPSACKS);

	for(int k=0;k<MULTI_KP_instance->N_KNAPSACKS;k++)
	{
		MULTI_KP_instance->MATRIX_WEIGHTS[k]=(int*)malloc(sizeof(int)*MULTI_KP_instance->N_ITEMS);
	}

	long int *sum_weight=new long int[MULTI_KP_instance->N_KNAPSACKS];


	if(MULTI_KP_instance->category==1)
	{

		cout << "Uncorrelated\n";

		for(int i=0;i<MULTI_KP_instance->N_ITEMS;i++)
		{
			MULTI_KP_instance->MULTI_PROFITS[i]=randomBETWEEN(1,MULTI_KP_instance->R);
		}

		for(int j=0;j<MULTI_KP_instance->N_KNAPSACKS;j++)
		{
			sum_weight[j]=0;

			for(int i=0;i<MULTI_KP_instance->N_ITEMS;i++)
			{

				MULTI_KP_instance->MATRIX_WEIGHTS[j][i]=randomBETWEEN(1,MULTI_KP_instance->R);

				sum_weight[j]+=MULTI_KP_instance->MATRIX_WEIGHTS[j][i];
			}
		}
	}


	if(MULTI_KP_instance->category==2)
	{

		cout <<  "\nWeakly correlated (weight and profit inverted respected to pisinger)\n";

		for(int i=0;i<MULTI_KP_instance->N_ITEMS;i++)
		{
			MULTI_KP_instance->MULTI_PROFITS[i]=randomBETWEEN(1,MULTI_KP_instance->R);
		}

		for(int j=0;j<MULTI_KP_instance->N_KNAPSACKS;j++)
		{
			sum_weight[j]=0;

			for(int i=0;i<MULTI_KP_instance->N_ITEMS;i++)
			{

				int LB_P=max( 1.0 , MULTI_KP_instance->MULTI_PROFITS[i] - (MULTI_KP_instance->R / 10.0) );

				int UB_P= MULTI_KP_instance->MULTI_PROFITS[i] + (MULTI_KP_instance->R / 10.0);

				MULTI_KP_instance->MATRIX_WEIGHTS[j][i]=randomBETWEEN(LB_P,UB_P);

				sum_weight[j]+=MULTI_KP_instance->MATRIX_WEIGHTS[j][i];
			}
		}
	}


	if(MULTI_KP_instance->category==3)
	{

		cout <<  "\nStrongly correlated\n";

		cout << "NOT IMPLEMENTED\n";

		exit(-1);

	}

	if(MULTI_KP_instance->category==4)
	{
		cout <<  "\nInverse strongly correlated\n";

		cout << "NOT IMPLEMENTED\n";

		exit(-1);

	}


	if(MULTI_KP_instance->category==5)
	{

		cout <<  "\nAlmost strongly correlated (weight and profit inverted respected to pisinger)\n";

		for(int i=0;i<MULTI_KP_instance->N_ITEMS;i++)
		{
			MULTI_KP_instance->MULTI_PROFITS[i]=randomBETWEEN(1,MULTI_KP_instance->R);
		}

		for(int j=0;j<MULTI_KP_instance->N_KNAPSACKS;j++)
		{
			sum_weight[j]=0;

			for(int i=0;i<MULTI_KP_instance->N_ITEMS;i++)
			{

				int LB_P= MULTI_KP_instance->MULTI_PROFITS[i] + (MULTI_KP_instance->R / 10.0) - (MULTI_KP_instance->R / 500.0);

				int UB_P= MULTI_KP_instance->MULTI_PROFITS[i] + (MULTI_KP_instance->R / 10.0) + (MULTI_KP_instance->R / 500.0);

				MULTI_KP_instance->MATRIX_WEIGHTS[j][i]=randomBETWEEN(LB_P,UB_P);

				sum_weight[j]+=MULTI_KP_instance->MATRIX_WEIGHTS[j][i];
			}
		}
	}


	if(MULTI_KP_instance->category==6)
	{

		cout << "Subset-sum\n";

		if(MULTI_KP_instance->N_KNAPSACKS!=1)
		{
			cout << "\n\n***NOT POSSIBLE WITH MULTIPLE KP***\n\n";
			exit(-1);
		}

		for(int j=0;j<MULTI_KP_instance->N_KNAPSACKS;j++)
		{
			sum_weight[j]=0;

			for(int i=0;i<MULTI_KP_instance->N_ITEMS;i++)
			{


				int dummy=randomBETWEEN(1,MULTI_KP_instance->R);

				MULTI_KP_instance->MULTI_PROFITS[i]=dummy;

				MULTI_KP_instance->MATRIX_WEIGHTS[j][i]=dummy;

				sum_weight[j]+=MULTI_KP_instance->MATRIX_WEIGHTS[j][i];
			}
		}
	}

	MULTI_KP_instance->RHS_CAPACITIES=(int*)malloc(sizeof(int)*MULTI_KP_instance->N_KNAPSACKS);

	for(int j=0;j<MULTI_KP_instance->N_KNAPSACKS;j++)
	{
		MULTI_KP_instance->RHS_CAPACITIES[j]= (int) (sum_weight[j] * MULTI_KP_instance->perc_cap/100.0);;
	}

	delete []sum_weight;

	for(int j=0;j<MULTI_KP_instance->N_KNAPSACKS;j++)
	{

		for(int i=0;i<MULTI_KP_instance->N_ITEMS;i++)
		{
			if(MULTI_KP_instance->MATRIX_WEIGHTS[j][i]>= MULTI_KP_instance->RHS_CAPACITIES[j])
			{
				MULTI_KP_instance->MATRIX_WEIGHTS[j][i]= MULTI_KP_instance->RHS_CAPACITIES[j];
			}
		}
	}


	////////////////////////////////////////////////////////////////////////////////////////////////////
	MULTI_KP_instance->n_conflicts=0;

	MULTI_KP_instance->CONF_MATRIX=new int*[MULTI_KP_instance->N_ITEMS];
	for ( int j = 0; j < MULTI_KP_instance->N_ITEMS; j++ )
	{
		MULTI_KP_instance->CONF_MATRIX[j]=new int[MULTI_KP_instance->N_ITEMS];
	}

	for ( int j = 0; j < MULTI_KP_instance->N_ITEMS; j++ )
	{
		for ( int i = 0; i < MULTI_KP_instance->N_ITEMS; i++ )
		{
			MULTI_KP_instance->CONF_MATRIX[j][i]=0;
		}
	}

	for(int i=0; i<MULTI_KP_instance->N_ITEMS; i++)
	{

		for(int j=i+1; j<MULTI_KP_instance->N_ITEMS; j++)
		{

			double n_random=randomBETWEEN_double(0,1);

			if(n_random<=MULTI_KP_instance->conflict_density)
			{

				MULTI_KP_instance->CONF_MATRIX[i][j]=1;
				MULTI_KP_instance->CONF_MATRIX[j][i]=1;

				MULTI_KP_instance->n_conflicts++;
			}
		}
	}
	////////////////////////////////////////////////////////////////////////////////////////////////////


#ifdef print_ist_features


	cout << "MATRIX_WEIGHTS\n";
	for(int j=0;j<MULTI_KP_instance->N_KNAPSACKS;j++)
	{
		for(int i=0;i<MULTI_KP_instance->N_ITEMS;i++)
		{
			cout <<  MULTI_KP_instance->MATRIX_WEIGHTS[j][i] << "\t";
		}
		cout << endl;
	}
	cout << endl;

	cout << "RHS_CAPACITIES\n";
	for(int j=0;j<MULTI_KP_instance->N_KNAPSACKS;j++)
	{
		cout <<  MULTI_KP_instance->RHS_CAPACITIES[j] << "\t";
	}
	cout << endl;

	cout << "MULTI_PROFITS\n";
	for(int i=0;i<MULTI_KP_instance->N_ITEMS;i++)
	{
		cout <<  MULTI_KP_instance->MULTI_PROFITS[i] << "\t";
	}
	cout << endl;

#endif


	cout << "Instance generated!\n";
	///////////////////////////////////////////////////////////

}

/*********************************/
void free_data_prob_multi_kp(data_MULTI_KP *MULTI_KP_instance)
/*********************************/
{

	cout << "FREEEEEE\n\n";

	for ( int j = 0; j < MULTI_KP_instance->N_ITEMS; j++ )
	{
		delete []MULTI_KP_instance->CONF_MATRIX[j];
	}
	delete []MULTI_KP_instance->CONF_MATRIX;


	free(MULTI_KP_instance->MULTI_PROFITS);

	for(int k=0;k<MULTI_KP_instance->N_KNAPSACKS;k++)
	{
		free(MULTI_KP_instance->MATRIX_WEIGHTS[k]);
	}
	free(MULTI_KP_instance->MATRIX_WEIGHTS);


	free(MULTI_KP_instance->RHS_CAPACITIES);

}


/***********************************************************************************/
void multi_kp_load_cplex(data_MULTI_KP *MULTI_KP_instance)
/***********************************************************************************/
{

	int i,j;


	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// * setting the CPLEX environment

	//opening the environment
	MULTI_KP_instance->env_MULTI_KP=CPXopenCPLEX(&(MULTI_KP_instance->status));
	if(MULTI_KP_instance->status!=0)
	{
		printf("cannot open CPLEX environment\n");
		exit(-1);
	}

	//opening the pointer to the problem
	MULTI_KP_instance->lp_MULTI_KP=CPXcreateprob(MULTI_KP_instance->env_MULTI_KP,&(MULTI_KP_instance->status),"KP");
	if(MULTI_KP_instance->status!=0)
	{
		printf("cannot create problem\n");
		exit(-1);
	}
	/////////////////////////////////////////////////////////////////////////////////////////////////////////

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// * creating the variables *
	MULTI_KP_instance->ccnt=MULTI_KP_instance->N_ITEMS;
	MULTI_KP_instance->obj=(double*) calloc(MULTI_KP_instance->ccnt,sizeof(double));
	MULTI_KP_instance->lb=(double*) calloc(MULTI_KP_instance->ccnt,sizeof(double));
	MULTI_KP_instance->ub=(double*) calloc(MULTI_KP_instance->ccnt,sizeof(double));
	MULTI_KP_instance->xctype=(char*) calloc(MULTI_KP_instance->ccnt,sizeof(char));


	MULTI_KP_instance->colname=(char**) calloc(MULTI_KP_instance->ccnt,sizeof(char*));
	for(i=0;i<MULTI_KP_instance->ccnt;i++){MULTI_KP_instance->colname[i]=(char*) calloc(2000,sizeof(char));}

	for(i=0; i<MULTI_KP_instance->ccnt; i++)
	{
		MULTI_KP_instance->obj[i]=MULTI_KP_instance->MULTI_PROFITS[i];
		MULTI_KP_instance->lb[i]=0.0;
		MULTI_KP_instance->ub[i]=1.0;
		MULTI_KP_instance->xctype[i]='B';
		sprintf(MULTI_KP_instance->colname[i], "x%d",i+1);
	}

	MULTI_KP_instance->status=CPXnewcols(MULTI_KP_instance->env_MULTI_KP,MULTI_KP_instance->lp_MULTI_KP,MULTI_KP_instance->ccnt,MULTI_KP_instance->obj,MULTI_KP_instance->lb,MULTI_KP_instance->ub,MULTI_KP_instance->xctype,MULTI_KP_instance->colname);
	if(MULTI_KP_instance->status!=0)
	{
		printf("error in CPXnewcols\n");
		exit(-1);
	}

	free(MULTI_KP_instance->obj);
	free(MULTI_KP_instance->lb);
	free(MULTI_KP_instance->ub);
	free(MULTI_KP_instance->xctype);

	for(i=0;i<MULTI_KP_instance->ccnt;i++){free(MULTI_KP_instance->colname[i]);}
	free(MULTI_KP_instance->colname);


	// * setting the objective function in the maximization form
	CPXchgobjsen(MULTI_KP_instance->env_MULTI_KP,MULTI_KP_instance->lp_MULTI_KP,CPX_MAX);
	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	//
	//	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	// * creating the knapsack constraints *

	for(int j=0;j<MULTI_KP_instance->N_KNAPSACKS;j++)
	{
		MULTI_KP_instance->rcnt=1;
		MULTI_KP_instance->nzcnt=MULTI_KP_instance->N_ITEMS;
		MULTI_KP_instance->rhs=(double*) calloc(MULTI_KP_instance->rcnt,sizeof(double));
		MULTI_KP_instance->sense=(char*) calloc(MULTI_KP_instance->rcnt,sizeof(double));

		MULTI_KP_instance->rhs[0]=MULTI_KP_instance->RHS_CAPACITIES[j];
		MULTI_KP_instance->sense[0]='L';


		MULTI_KP_instance->rmatbeg=(int*) calloc(MULTI_KP_instance->rcnt,sizeof(int));
		MULTI_KP_instance->rmatind=(int*) calloc(MULTI_KP_instance->nzcnt,sizeof(int));
		MULTI_KP_instance->rmatval=(double*) calloc(MULTI_KP_instance->nzcnt,sizeof(double));

		for(i=0; i<MULTI_KP_instance->N_ITEMS; i++)
		{
			MULTI_KP_instance->rmatval[i]=MULTI_KP_instance->MATRIX_WEIGHTS[j][i];
			MULTI_KP_instance->rmatind[i]=i;
		}

		MULTI_KP_instance->rmatbeg[0]=0;

		MULTI_KP_instance->status=CPXaddrows(MULTI_KP_instance->env_MULTI_KP,MULTI_KP_instance->lp_MULTI_KP,0,MULTI_KP_instance->rcnt,MULTI_KP_instance->nzcnt,MULTI_KP_instance->rhs,MULTI_KP_instance->sense,MULTI_KP_instance->rmatbeg,MULTI_KP_instance->rmatind,MULTI_KP_instance->rmatval,NULL,NULL);
		if(MULTI_KP_instance->status!=0)
		{
			printf("error in CPXaddrows\n");
			exit(-1);
		}

		free(MULTI_KP_instance->rmatbeg);
		free(MULTI_KP_instance->rmatval);
		free(MULTI_KP_instance->rmatind);
		free(MULTI_KP_instance->rhs);
		free(MULTI_KP_instance->sense);
	}
	//	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	//

	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	for(i=0; i<MULTI_KP_instance->N_ITEMS; i++)
	{
		for(j=i+1; j<MULTI_KP_instance->N_ITEMS; j++)
		{
			if(MULTI_KP_instance->CONF_MATRIX[i][j]==1)
			{

				//cout << "item conflict\t" << i << "\t" << j << endl;

				MULTI_KP_instance->rcnt=1;
				MULTI_KP_instance->nzcnt=2;
				MULTI_KP_instance->rhs=(double*) calloc(MULTI_KP_instance->rcnt,sizeof(double));
				MULTI_KP_instance->sense=(char*) calloc(MULTI_KP_instance->rcnt,sizeof(double));

				MULTI_KP_instance->rhs[0]=1.0;
				MULTI_KP_instance->sense[0]='L';


				MULTI_KP_instance->rmatbeg=(int*) calloc(MULTI_KP_instance->rcnt,sizeof(int));
				MULTI_KP_instance->rmatind=(int*) calloc(MULTI_KP_instance->nzcnt,sizeof(int));
				MULTI_KP_instance->rmatval=(double*) calloc(MULTI_KP_instance->nzcnt,sizeof(double));

				MULTI_KP_instance->rmatval[0]=1.0;
				MULTI_KP_instance->rmatind[0]=i;

				MULTI_KP_instance->rmatval[1]=1.0;
				MULTI_KP_instance->rmatind[1]=j;

				MULTI_KP_instance->rmatbeg[0]=0;

				MULTI_KP_instance->status=CPXaddrows(MULTI_KP_instance->env_MULTI_KP,MULTI_KP_instance->lp_MULTI_KP,0,MULTI_KP_instance->rcnt,MULTI_KP_instance->nzcnt,MULTI_KP_instance->rhs,MULTI_KP_instance->sense,MULTI_KP_instance->rmatbeg,MULTI_KP_instance->rmatind,MULTI_KP_instance->rmatval,NULL,NULL);
				if(MULTI_KP_instance->status!=0)
				{
					printf("error in CPXaddrows\n");
					exit(-1);
				}

				free(MULTI_KP_instance->rmatbeg);
				free(MULTI_KP_instance->rmatval);
				free(MULTI_KP_instance->rmatind);
				free(MULTI_KP_instance->rhs);
				free(MULTI_KP_instance->sense);

			}
		}
	}
	////////////////////////////////////////////////////////////////////////////////////////////////////////////




	//
	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// * writing the created ILP model on a file *
	MULTI_KP_instance->status=CPXwriteprob(MULTI_KP_instance->env_MULTI_KP,MULTI_KP_instance->lp_MULTI_KP,"MULTI_kp.lp",NULL);
	if(MULTI_KP_instance->status!=0) {
		printf("error in CPXwriteprob\n");
		exit(-1);
	}
	/////////////////////////////////////////////////////////////////////////////////////////////////////////




}

/***********************************************************************************/
void multi_kp_free_cplex(data_MULTI_KP *MULTI_KP_instance)
/***********************************************************************************/
{


	MULTI_KP_instance->status=CPXfreeprob(MULTI_KP_instance->env_MULTI_KP,&(MULTI_KP_instance->lp_MULTI_KP));
	if(MULTI_KP_instance->status!=0) {printf("error in CPXfreeprob\n");exit(-1);}

	MULTI_KP_instance->status=CPXcloseCPLEX(&(MULTI_KP_instance->env_MULTI_KP));
	if(MULTI_KP_instance->status!=0) {printf("error in CPXcloseCPLEX\n");exit(-1);}

}

/***********************************************************************************/
double multi_kp_solve_cplex(data_MULTI_KP *MULTI_KP_instance)
/***********************************************************************************/
{



	int i,j;

	// * Set video output *
	MULTI_KP_instance->status = CPXsetintparam (MULTI_KP_instance->env_MULTI_KP, CPX_PARAM_SCRIND, CPX_ON);
	if (MULTI_KP_instance->status)
	{
		printf ("error for CPX_PARAM_SCRIND\n");
	}


	//	// * Set  tolerance *
	//	MULTI_KP_instance->status = CPXsetdblparam (MULTI_KP_instance->env_MULTI_KP, CPX_PARAM_EPAGAP, 0.0);
	//	if (MULTI_KP_instance->status)
	//	{
	//		printf ("error for CPX_PARAM_EPAGAP\n");
	//	}
	//
	//	// * Set  tolerance *
	//	MULTI_KP_instance->status = CPXsetdblparam (MULTI_KP_instance->env_MULTI_KP, CPX_PARAM_EPGAP, 0.0);
	//	if (MULTI_KP_instance->status)
	//	{
	//		printf ("error for CPX_PARAM_EPGAP\n");
	//	}
	//
	//
	//	// * Set mip tolerances integrality *
	//	MULTI_KP_instance->status = CPXsetdblparam (MULTI_KP_instance->env_MULTI_KP, CPX_PARAM_EPINT, 0.0);
	//	if (MULTI_KP_instance->status)
	//	{
	//		printf ("error for CPX_PARAM_EPINTP\n");
	//	}
	//
	//	// * Set Feasibility tolerance *
	//	MULTI_KP_instance->status = CPXsetdblparam (MULTI_KP_instance->env_MULTI_KP, CPX_PARAM_EPRHS, 1e-9);
	//	if (MULTI_KP_instance->status)
	//	{
	//		printf ("error for CPX_PARAM_EPRHS\n");
	//	}
	//
	// * Set number of CPU*

	MULTI_KP_instance->status = CPXsetintparam (MULTI_KP_instance->env_MULTI_KP, CPX_PARAM_THREADS, 1);
	if (MULTI_KP_instance->status)
	{
		printf ("error for CPX_PARAM_THREADS\n");
	}

	// * Set time limit *

	MULTI_KP_instance->status = CPXsetdblparam (MULTI_KP_instance->env_MULTI_KP, CPX_PARAM_TILIM,MULTI_KP_instance->timeLimit);
	if (MULTI_KP_instance->status)
	{
		printf ("error for CPX_PARAM_TILIM\n");
	}



	///////////////////////////////////////////////////////////////////////////////////
	// * solving the MIP model
	clock_t time_start=clock();


	MULTI_KP_instance->status=CPXmipopt(MULTI_KP_instance->env_MULTI_KP,MULTI_KP_instance->lp_MULTI_KP);
	if(MULTI_KP_instance->status!=0)
	{
		printf("error in CPXmipopt\n");
		exit(-1);
	}

	clock_t time_end=clock();
	MULTI_KP_instance->solution_time=(double)(time_end-time_start)/(double)CLOCKS_PER_SEC;

	printf("\n\nsolution_time ->\t\%f",MULTI_KP_instance->solution_time);


	///////////////////////////////////////////////////////////////////////////////////


	// * getting the solution

	MULTI_KP_instance->x=(double*) calloc(MULTI_KP_instance->N_ITEMS,sizeof(double));


	MULTI_KP_instance->status=CPXgetmipx(MULTI_KP_instance->env_MULTI_KP,MULTI_KP_instance->lp_MULTI_KP,MULTI_KP_instance->x,0,MULTI_KP_instance->N_ITEMS-1);
	if(MULTI_KP_instance->status!=0)
	{
		printf("error in CPXgetmipx\n");
		exit(-1);
	}

	MULTI_KP_instance->status=CPXgetmipobjval(MULTI_KP_instance->env_MULTI_KP,MULTI_KP_instance->lp_MULTI_KP,&(MULTI_KP_instance->objval));
	if(MULTI_KP_instance->status!=0)
	{
		printf("error in CPXgetmipobjval\n");
		exit(-1);
	}

	printf("\n\nMIP solution value ->\t\%f\n",MULTI_KP_instance->objval);

#ifdef print_sol
	printf("\n\nSolution\n");
	for (i = 0; i < MULTI_KP_instance->N_ITEMS; i++)
	{
		printf("item %d -> %d\n",i+1 ,(int)(MULTI_KP_instance->x[i]+0.5));
	}
	printf("\n\n");
#endif

	MULTI_KP_instance->status=CPXgetbestobjval(MULTI_KP_instance->env_MULTI_KP,MULTI_KP_instance->lp_MULTI_KP,&(MULTI_KP_instance->bestobjval));
	if(MULTI_KP_instance->status!=0)
	{
		printf("error in CPXgetbestobjval\n");
		exit(-1);
	}

	MULTI_KP_instance->stat=CPXgetstat(MULTI_KP_instance->env_MULTI_KP,MULTI_KP_instance->lp_MULTI_KP);
	MULTI_KP_instance->nodecount = CPXgetnodecnt(MULTI_KP_instance->env_MULTI_KP, MULTI_KP_instance->lp_MULTI_KP);


	//	///////////////////////////////////////////////////////////////////////////////////
	//	/* linear programming relaxation*/
	//
	MULTI_KP_instance->solution_time_lp=0;

	MULTI_KP_instance->status = CPXchgprobtype (MULTI_KP_instance->env_MULTI_KP, MULTI_KP_instance->lp_MULTI_KP, CPXPROB_LP);

	clock_t time_start_lp=clock();
	MULTI_KP_instance->status=CPXlpopt(MULTI_KP_instance->env_MULTI_KP,MULTI_KP_instance->lp_MULTI_KP);
	if(MULTI_KP_instance->status!=0)
	{
		printf("err_FILEor in CPXlpopt slave solve\n");
		exit(-1);
	}

	clock_t time_end_lp=clock();

	MULTI_KP_instance->solution_time_lp=(double)(time_end_lp-time_start_lp)/(double)CLOCKS_PER_SEC;

	MULTI_KP_instance->cplex_lp;

	MULTI_KP_instance->status=CPXgetobjval(MULTI_KP_instance->env_MULTI_KP,MULTI_KP_instance->lp_MULTI_KP,&MULTI_KP_instance->cplex_lp);
	if(MULTI_KP_instance->status!=0)
	{
		printf("error in CPXgetmipobjval\n");
		exit(-1);
	}

	printf("\n\nLP solution value ->\t\%f",MULTI_KP_instance->cplex_lp);

	MULTI_KP_instance->status=CPXgetx(MULTI_KP_instance->env_MULTI_KP,MULTI_KP_instance->lp_MULTI_KP,MULTI_KP_instance->x,0,MULTI_KP_instance->N_ITEMS-1);
	if(MULTI_KP_instance->status!=0)
	{
		printf("error in CPXgetmipx\n");
		exit(-1);
	}

#ifdef print_sol_frac
	printf("\n\nSolution\n");
	for (i = 0; i < MULTI_KP_instance->N_ITEMS; i++)
	{
		printf("item %d -> %.3f\n",i+1,MULTI_KP_instance->x[i]);
	}
	printf("\n");
#endif

	///////////////////////////////////////////////////////////////////////////////////
	//
	MULTI_KP_instance->cur_numcols=CPXgetnumcols(MULTI_KP_instance->env_MULTI_KP,MULTI_KP_instance->lp_MULTI_KP);
	MULTI_KP_instance->cur_numrows=CPXgetnumrows(MULTI_KP_instance->env_MULTI_KP,MULTI_KP_instance->lp_MULTI_KP);


	cout << "\n\nDONE\n";

	free(MULTI_KP_instance->x);




	return MULTI_KP_instance->objval;

}


/***********************************************************************************/
void add_cover_extended(data_MULTI_KP *MULTI_KP_instance,double *cover, double RHS)
/***********************************************************************************/
{

	int i,j;

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// * creating the knapsack constraint *
	MULTI_KP_instance->rcnt=1;
	MULTI_KP_instance->nzcnt=MULTI_KP_instance->N_ITEMS;
	MULTI_KP_instance->rhs=(double*) calloc(MULTI_KP_instance->rcnt,sizeof(double));
	MULTI_KP_instance->sense=(char*) calloc(MULTI_KP_instance->rcnt,sizeof(double));

	MULTI_KP_instance->sense[0]='L';


	MULTI_KP_instance->rmatbeg=(int*) calloc(MULTI_KP_instance->rcnt,sizeof(int));
	MULTI_KP_instance->rmatind=(int*) calloc(MULTI_KP_instance->nzcnt,sizeof(int));
	MULTI_KP_instance->rmatval=(double*) calloc(MULTI_KP_instance->nzcnt,sizeof(double));

	int counter=0;
	for(i=0; i<MULTI_KP_instance->N_ITEMS; i++)
	{
		if(cover[i]>0.5)
		{
			MULTI_KP_instance->rmatval[i]=1.0;
			counter++;
		}
		else{
			MULTI_KP_instance->rmatval[i]=0.0;
		}
		MULTI_KP_instance->rmatind[i]=i;
	}

	MULTI_KP_instance->rhs[0]=RHS-1;


	MULTI_KP_instance->rmatbeg[0]=0;

	MULTI_KP_instance->status=CPXaddrows(MULTI_KP_instance->env_MULTI_KP,MULTI_KP_instance->lp_MULTI_KP,0,MULTI_KP_instance->rcnt,MULTI_KP_instance->nzcnt,MULTI_KP_instance->rhs,MULTI_KP_instance->sense,MULTI_KP_instance->rmatbeg,MULTI_KP_instance->rmatind,MULTI_KP_instance->rmatval,NULL,NULL);
	if(MULTI_KP_instance->status!=0)
	{
		printf("error in CPXaddrows\n");
		exit(-1);
	}

	free(MULTI_KP_instance->rmatbeg);
	free(MULTI_KP_instance->rmatval);
	free(MULTI_KP_instance->rmatind);
	free(MULTI_KP_instance->rhs);
	free(MULTI_KP_instance->sense);
	/////////////////////////////////////////////////////////////////////////////////////////////////////////


	//		/////////////////////////////////////////////////////////////////////////////////////////////////////////
	//		// * writing the created ILP model on a file *
	//		MULTI_KP_instance->status=CPXwriteprob(MULTI_KP_instance->env_MULTI_KP,MULTI_KP_instance->lp_MULTI_KP,"sep_cover.lp",NULL);
	//		if(MULTI_KP_instance->status!=0) {
	//			printf("error in CPXwriteprob\n");
	//			exit(-1);
	//		}
	//		exit(-1);
	//		/////////////////////////////////////////////////////////////////////////////////////////////////////////

}


/***********************************************************************************/
void add_cover(data_MULTI_KP *MULTI_KP_instance,double *cover)
/***********************************************************************************/
{

	int i,j;

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// * creating the knapsack constraint *
	MULTI_KP_instance->rcnt=1;
	MULTI_KP_instance->nzcnt=MULTI_KP_instance->N_ITEMS;
	MULTI_KP_instance->rhs=(double*) calloc(MULTI_KP_instance->rcnt,sizeof(double));
	MULTI_KP_instance->sense=(char*) calloc(MULTI_KP_instance->rcnt,sizeof(double));

	MULTI_KP_instance->sense[0]='L';


	MULTI_KP_instance->rmatbeg=(int*) calloc(MULTI_KP_instance->rcnt,sizeof(int));
	MULTI_KP_instance->rmatind=(int*) calloc(MULTI_KP_instance->nzcnt,sizeof(int));
	MULTI_KP_instance->rmatval=(double*) calloc(MULTI_KP_instance->nzcnt,sizeof(double));

	int counter=0;
	for(i=0; i<MULTI_KP_instance->N_ITEMS; i++)
	{
		if(cover[i]>0.5)
		{
			MULTI_KP_instance->rmatval[i]=1.0;
			counter++;
		}
		else{
			MULTI_KP_instance->rmatval[i]=0.0;
		}
		MULTI_KP_instance->rmatind[i]=i;
	}

	MULTI_KP_instance->rhs[0]=counter-1;


	MULTI_KP_instance->rmatbeg[0]=0;

	MULTI_KP_instance->status=CPXaddrows(MULTI_KP_instance->env_MULTI_KP,MULTI_KP_instance->lp_MULTI_KP,0,MULTI_KP_instance->rcnt,MULTI_KP_instance->nzcnt,MULTI_KP_instance->rhs,MULTI_KP_instance->sense,MULTI_KP_instance->rmatbeg,MULTI_KP_instance->rmatind,MULTI_KP_instance->rmatval,NULL,NULL);
	if(MULTI_KP_instance->status!=0)
	{
		printf("error in CPXaddrows\n");
		exit(-1);
	}

	free(MULTI_KP_instance->rmatbeg);
	free(MULTI_KP_instance->rmatval);
	free(MULTI_KP_instance->rmatind);
	free(MULTI_KP_instance->rhs);
	free(MULTI_KP_instance->sense);
	/////////////////////////////////////////////////////////////////////////////////////////////////////////


	//	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	// * writing the created ILP model on a file *
	//	MULTI_KP_instance->status=CPXwriteprob(MULTI_KP_instance->env_MULTI_KP,MULTI_KP_instance->lp_MULTI_KP,"sep_cover.lp",NULL);
	//	if(MULTI_KP_instance->status!=0) {
	//		printf("error in CPXwriteprob\n");
	//		exit(-1);
	//	}
	//	/////////////////////////////////////////////////////////////////////////////////////////////////////////

}


/***********************************************************************************/
double multi_kp_solve_cplex_LP(data_MULTI_KP *MULTI_KP_instance,double *point)
/***********************************************************************************/
{

	int i,j;

	// * Set video output *
	MULTI_KP_instance->status = CPXsetintparam (MULTI_KP_instance->env_MULTI_KP, CPX_PARAM_SCRIND, CPX_OFF);
	if (MULTI_KP_instance->status)
	{
		printf ("error for CPX_PARAM_SCRIND\n");
	}

	// * Set number of CPU*

	MULTI_KP_instance->status = CPXsetintparam (MULTI_KP_instance->env_MULTI_KP, CPX_PARAM_THREADS, 1);
	if (MULTI_KP_instance->status)
	{
		printf ("error for CPX_PARAM_THREADS\n");
	}

	//	// * Set time limit *
	//
	//	KP_instance->status = CPXsetdblparam (KP_instance->env, CPX_PARAM_TILIM,KP_instance->timeLimit);
	//	if (KP_instance->status)
	//	{
	//		printf ("error for CPX_PARAM_TILIM\n");
	//	}


	///////////////////////////////////////////////////////////////////////////////////
	/* linear programming relaxation*/

	double solution_time_lp=0;

	MULTI_KP_instance->status = CPXchgprobtype (MULTI_KP_instance->env_MULTI_KP, MULTI_KP_instance->lp_MULTI_KP, CPXPROB_LP);

	clock_t time_start_lp=clock();
	MULTI_KP_instance->status=CPXlpopt(MULTI_KP_instance->env_MULTI_KP,MULTI_KP_instance->lp_MULTI_KP);
	if(MULTI_KP_instance->status!=0)
	{
		printf("err_FILEor in CPXlpopt slave solve\n");
		exit(-1);
	}

	clock_t time_end_lp=clock();

	solution_time_lp=(double)(time_end_lp-time_start_lp)/(double)CLOCKS_PER_SEC;

	double cplex_lp;

	MULTI_KP_instance->status=CPXgetobjval(MULTI_KP_instance->env_MULTI_KP,MULTI_KP_instance->lp_MULTI_KP,&cplex_lp);
	if(MULTI_KP_instance->status!=0)
	{
		printf("error in CPXgetmipobjval\n");
		exit(-1);
	}

	printf("LP solution value ->\t\%f\t",cplex_lp);

	MULTI_KP_instance->status=CPXgetx(MULTI_KP_instance->env_MULTI_KP,MULTI_KP_instance->lp_MULTI_KP,point,0,MULTI_KP_instance->N_ITEMS-1);
	if(MULTI_KP_instance->status!=0)
	{
		printf("error in CPXgetmipx\n");
		exit(-1);
	}

	//		printf("\n\nSolution\n");
	//		for (i = 0; i <MULTI_KP_instance->N_ITEMS; i++)
	//		{
	//			printf("item %d -> %.3f \n",i+1,point[i]);
	//		}
	//		printf("\n");
	//		cin.get();

	///////////////////////////////////////////////////////////////////////////////////

	//	KP_instance->cur_numcols=CPXgetnumcols(KP_instance->env,KP_instance->lp);
	//	KP_instance->cur_numrows=CPXgetnumrows(KP_instance->env,KP_instance->lp);
	//
	//	printf("\nnumcols\t%d\n",KP_instance->cur_numcols);
	//	printf("\nnumrows\t%d\n",KP_instance->cur_numrows);



	return cplex_lp;

}





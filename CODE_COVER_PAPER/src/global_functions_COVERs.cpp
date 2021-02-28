

#include "global_functions_COVERs.h"

/***********************************************************************************/
int compare(const void *p, const void *q)
/***********************************************************************************/
{
	//SORT BY NON INCREASING SCORES
	double l = ((struct valuesSTR *)p)->score;
	double r = ((struct valuesSTR *)q)->score;
	return ((l < r)  ? 1:-1);
}


/***********************************************************************************/
int extend_cover(int cover_size,int cap,double *cover,int *weight)
/***********************************************************************************/
{

	//		for (int i = 0; i < cover_size; i++)
	//		{
	//			cout << cover[i] << "\t weight \t" <<weight[i] << endl;
	//		}
	//		cout << endl;

	int lifted_items=0;

	int max_weight=0;
	for (int i = 0; i < cover_size; i++)
	{
		if(max_weight<weight[i] && cover[i]>0.5)
		{
			max_weight=weight[i];
		}
	}

	//	cout << "max_weight\t" << max_weight << endl;

	for (int i = 0; i < cover_size; i++)
	{
		if(cover[i]<0.5 && weight[i] >= max_weight)
		{
			cover[i]=1.0;
			lifted_items++;

			//	cout << "lifted->\t" << i << endl;
		}
	}

	//	cin.get();

	return lifted_items;

}

/***********************************************************************************/
double solve_DP_lex(data_COVER *COVER_instance,int index)
/***********************************************************************************/
{

	double EPSILON_MAGIC= 1E-4 / (3*COVER_instance->cover_size);

	//EPSILON_MAGIC=0;

	int i,j;

	for (i = 0; i < COVER_instance->cover_size; i++)
	{
		COVER_instance->profit_sep_magic[i]= COVER_instance->profit_sep[i] + EPSILON_MAGIC ;
	}


	///////////////////////////////////////////////////////////////////////////////////
	// * solving the DP
	clock_t time_start=clock();


	COVER_instance->objval_cover=DP_kp01_advanced(COVER_instance->cover_size, COVER_instance->cap_sep[index] , COVER_instance->profit_sep_magic,COVER_instance->weight_sep[index],COVER_instance->cover, true , COVER_instance->order_DP[index]);

	clock_t time_end=clock();
	double solution_time=(double)(time_end-time_start)/(double)CLOCKS_PER_SEC;

	//printf("\n\nsolution_time ->\t\%f",solution_time);


	///////////////////////////////////////////////////////////////////////////////////

	double opt_val =   COVER_instance->cover_size;

	for (i = 0; i < COVER_instance->cover_size; i++)
	{

		opt_val -= COVER_instance->point[i];

		opt_val -= COVER_instance->cover[i] * COVER_instance->profit_sep[i];

		COVER_instance->cover[i]=1-COVER_instance->cover[i];
	}


	return opt_val;

}

/******************************************/
int make_maximal(data_COVER *COVER_instance,int cover_size,int cap_sep,double *cover,int *weight_sep,double *profit_sep,int index)
/******************************************/
{

	int items_added=0;

	double available_cap=cap_sep;

	for (int i = 0; i < cover_size; i++)
	{
		available_cap-=(cover[i]*weight_sep[i]);
	}

	//cout << "\navailable_cap\t" << available_cap << endl;


	for (int ii = 0; ii < cover_size; ii++)
	{

		int i=COVER_instance->order_MAXIMALITY[index][ii];

		if(cover[i]>0){continue;}

		if(weight_sep[i]<=available_cap)
		{

			//	cout << "item\t" << i << "\t anti cover \t"<< cover[i] << "\t p \t" << profit_sep[i] << "\t w \t" << weight_sep[i] << "\t cap avail.\t " << available_cap <<endl;

			available_cap-=weight_sep[i];

			cover[i]=1;

			items_added++;

			if(profit_sep[i]>0.000001)
			{
				cout << "profit_sep item\t" << i << "\t" << profit_sep[i] << "\tnot maximal\n";
			}
		}
	}

	//cout << "items_added\t" << items_added<< endl;


	return items_added;

}


/***********************************************************************************/
int compute_max_i(int *weights,int size,int cap_sep)
/***********************************************************************************/
{
	int max_i=0;

	vector < int > weights_sorted;

	for(int j=0; j< size; j++)
	{
		weights_sorted.push_back(weights[j]);
	}

	sort(weights_sorted.begin(), weights_sorted.end());

	int dummy_weight=0;

	for(int j=0; j<size; j++)
	{

		dummy_weight+=weights_sorted[j];

		if(dummy_weight>cap_sep)
		{
			break;
		}

		max_i=j+1;
	}

	return max_i;
}

/***********************************************************************************/
void  COVER_DATA_init(data_COVER *COVER_instance,int cover_size,double *profits,int **weights, int *capacity,int n_KP,int shuffle)
/***********************************************************************************/
{

	COVER_instance->cover_size=cover_size;


	COVER_instance->capacity=(int*) calloc(n_KP,sizeof(int));

	for(int i=0; i<n_KP; i++)
	{
		COVER_instance->capacity[i]=capacity[i];
	}

	COVER_instance->magic_sum=(long int*) calloc(COVER_instance->cover_size,sizeof(long int));
	COVER_instance->profit_sep_magic=(double*) calloc(COVER_instance->cover_size,sizeof(double));


	for(int i=0; i<COVER_instance->cover_size; i++)
	{
		COVER_instance->magic_sum[i]=0;

	}

	if(shuffle==1)
	{
		COVER_instance->order_DP.resize(n_KP);

		for (int j = 0; j < n_KP; j++)
		{
			COVER_instance->order_DP[j].resize(COVER_instance->cover_size);

			for (int jj = 0; jj < COVER_instance->cover_size; jj++)
			{
				COVER_instance->order_DP[j][jj]=jj;
			}
		}
	}

	if(shuffle==2)
	{

		cout << "\n->NON-increasing weight\n\n";

		COVER_instance->order_DP.resize(n_KP);

		for (int j = 0; j < n_KP; j++)
		{

			COVER_instance->order_DP[j].resize(COVER_instance->cover_size);

			valuesSTR *values_all=new valuesSTR[COVER_instance->cover_size];

			for (int jj = 0; jj < COVER_instance->cover_size; jj++)
			{

				values_all[jj].item=jj;
				values_all[jj].score=weights[j][jj];
			}

			qsort (values_all,COVER_instance->cover_size, sizeof(valuesSTR), compare);

			for (int jj = 0; jj < COVER_instance->cover_size; jj++)
			{

				//cout << "item\t" << values_all[jj].item << "\t" << values_all[jj].score << endl;

				COVER_instance->order_DP[j][jj]=values_all[jj].item;
			}
			//cout << endl;
			//cin.get();

			delete []values_all;


		}

	}


	if(shuffle==3)
	{

		cout << "\n->NON-decreasing weight\n\n";

		COVER_instance->order_DP.resize(n_KP);

		for (int j = 0; j < n_KP; j++)
		{

			COVER_instance->order_DP[j].resize(COVER_instance->cover_size);

			valuesSTR *values_all=new valuesSTR[COVER_instance->cover_size];

			for (int jj = 0; jj < COVER_instance->cover_size; jj++)
			{

				values_all[jj].item=jj;
				values_all[jj].score=1.0/weights[j][jj];
			}

			qsort (values_all,COVER_instance->cover_size, sizeof(valuesSTR), compare);

			for (int jj = 0; jj < COVER_instance->cover_size; jj++)
			{

				//cout << "item\t" << values_all[jj].item << "\t" << values_all[jj].score << "\t" << 1.0/values_all[jj].score << endl;

				COVER_instance->order_DP[j][jj]=values_all[jj].item;
			}
			//cout << endl;
			//cin.get();

			delete []values_all;
		}

	}

	if(shuffle==4)
	{

		cout << "\n->PROFIT OVER WEIGTH RATIO\n\n";

		COVER_instance->order_DP.resize(n_KP);

		for (int j = 0; j < n_KP; j++)
		{

			COVER_instance->order_DP[j].resize(COVER_instance->cover_size);

			valuesSTR *values_all=new valuesSTR[COVER_instance->cover_size];

			for (int jj = 0; jj < COVER_instance->cover_size; jj++)
			{

				values_all[jj].item=jj;
				values_all[jj].score=profits[jj]/weights[j][jj];
			}

			qsort (values_all,COVER_instance->cover_size, sizeof(valuesSTR), compare);

			for (int jj = 0; jj < COVER_instance->cover_size; jj++)
			{

				//cout << "item\t" << values_all[jj].item << "\t" << values_all[jj].score << endl;

				COVER_instance->order_DP[j][jj]=values_all[jj].item;
			}

			delete []values_all;
		}

	}


	//////////////////////////////////////////////////////////////////////////////////////////////////////
	cout << "\n->NON-decreasing weight FOR MAXIMALITY\n\n";

	COVER_instance->order_MAXIMALITY.resize(n_KP);

	for (int j = 0; j < n_KP; j++)
	{

		COVER_instance->order_MAXIMALITY[j].resize(COVER_instance->cover_size);

		valuesSTR *values_all=new valuesSTR[COVER_instance->cover_size];

		for (int jj = 0; jj < COVER_instance->cover_size; jj++)
		{

			values_all[jj].item=jj;
			values_all[jj].score=1.0/weights[j][jj];
		}

		qsort (values_all,COVER_instance->cover_size, sizeof(valuesSTR), compare);

		for (int jj = 0; jj < COVER_instance->cover_size; jj++)
		{

			//cout << "item\t" << values_all[jj].item << "\t" << values_all[jj].score << endl;

			COVER_instance->order_MAXIMALITY[j][jj]=values_all[jj].item;
		}
		//cout << endl;
		//cin.get();

		delete []values_all;
	}
	//////////////////////////////////////////////////////////////////////////////////////////////////////


	COVER_instance->point=(double*) calloc(COVER_instance->cover_size,sizeof(double));
	COVER_instance->profit_sep=(double*) calloc(COVER_instance->cover_size,sizeof(double));

	COVER_instance->weight_sep=(int**) calloc(COVER_instance->cover_size,sizeof(int*));
	for(int i=0; i<n_KP; i++)
	{
		COVER_instance->weight_sep[i]=(int*) calloc(COVER_instance->cover_size,sizeof(int));

	}

	COVER_instance->cover=(double*) calloc(COVER_instance->cover_size,sizeof(double));

	COVER_instance->cap_sep=(int*) calloc(n_KP,sizeof(int));
	COVER_instance->max_i=(int*) calloc(n_KP,sizeof(int));

	for(int i=0; i<n_KP; i++)
	{
		COVER_instance->cap_sep[i]=0;

		for(int j=0; j<COVER_instance->cover_size; j++)
		{
			COVER_instance->weight_sep[i][j]=weights[i][j];
			COVER_instance->cap_sep[i]+=weights[i][j];
		}
		COVER_instance->cap_sep[i]-=(capacity[i]+1);


		COVER_instance->max_i[i]=compute_max_i(weights[i],COVER_instance->cover_size,COVER_instance->cap_sep[i]);
	}

}


/***********************************************************************************/
void  COVER_DATA_free(data_COVER *COVER_instance,int n_KP)
/***********************************************************************************/
{


	free(COVER_instance->capacity);

	free(COVER_instance->magic_sum);
	free(COVER_instance->profit_sep_magic);
	free(COVER_instance->profit_sep);

	for(int i=0; i<n_KP; i++)
	{
		free(COVER_instance->weight_sep[i]);
	}
	free(COVER_instance->weight_sep);

	free(COVER_instance->cover);
	free(COVER_instance->point);

	free(COVER_instance->cap_sep);
	free(COVER_instance->max_i);

}

/***********************************************************************************/
void load_cplex_cover(data_COVER *COVER_instance,int index,int *weights, int capacity)
/***********************************************************************************/
{

	int i,j;


	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// * setting the CPLEX environment

	//opening the environment
	COVER_instance->env_sep[index]=CPXopenCPLEX(&(COVER_instance->status));
	if(COVER_instance->status!=0)
	{
		printf("cannot open CPLEX environment\n");
		exit(-1);
	}

	//opening the pointer to the problem
	COVER_instance->lp_sep[index]=CPXcreateprob(COVER_instance->env_sep[index],&(COVER_instance->status),"SEP");
	if(COVER_instance->status!=0)
	{
		printf("cannot create problem\n");
		exit(-1);
	}
	/////////////////////////////////////////////////////////////////////////////////////////////////////////

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// * creating the variables *
	COVER_instance->ccnt=COVER_instance->cover_size;
	COVER_instance->obj=(double*) calloc(COVER_instance->ccnt,sizeof(double));
	COVER_instance->lb=(double*) calloc(COVER_instance->ccnt,sizeof(double));
	COVER_instance->ub=(double*) calloc(COVER_instance->ccnt,sizeof(double));
	COVER_instance->xctype=(char*) calloc(COVER_instance->ccnt,sizeof(char));


	COVER_instance->colname=(char**) calloc(COVER_instance->ccnt,sizeof(char*));
	for(i=0;i<COVER_instance->ccnt;i++){COVER_instance->colname[i]=(char*) calloc(2000,sizeof(char));}

	for(i=0; i<COVER_instance->ccnt; i++)
	{
		COVER_instance->obj[i]=0.0;
		COVER_instance->lb[i]=0.0;
		COVER_instance->ub[i]=1.0;
		COVER_instance->xctype[i]='B';
		sprintf(COVER_instance->colname[i], "x%d",i+1);
	}

	COVER_instance->status=CPXnewcols(COVER_instance->env_sep[index],COVER_instance->lp_sep[index],COVER_instance->ccnt,COVER_instance->obj,COVER_instance->lb,COVER_instance->ub,COVER_instance->xctype,COVER_instance->colname);
	if(COVER_instance->status!=0)
	{
		printf("error in CPXnewcols\n");
		exit(-1);
	}

	free(COVER_instance->obj);
	free(COVER_instance->lb);
	free(COVER_instance->ub);
	free(COVER_instance->xctype);

	for(i=0;i<COVER_instance->ccnt;i++){free(COVER_instance->colname[i]);}
	free(COVER_instance->colname);


	// * setting the objective function in the maximization form
	CPXchgobjsen(COVER_instance->env_sep[index],COVER_instance->lp_sep[index],CPX_MAX);
	/////////////////////////////////////////////////////////////////////////////////////////////////////////

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// * creating the knapsack constraint *
	COVER_instance->rcnt=1;
	COVER_instance->nzcnt=COVER_instance->cover_size;
	COVER_instance->rhs=(double*) calloc(COVER_instance->rcnt,sizeof(double));
	COVER_instance->sense=(char*) calloc(COVER_instance->rcnt,sizeof(double));

	COVER_instance->rhs[0]=capacity+1;
	COVER_instance->sense[0]='G';


	COVER_instance->rmatbeg=(int*) calloc(COVER_instance->rcnt,sizeof(int));
	COVER_instance->rmatind=(int*) calloc(COVER_instance->nzcnt,sizeof(int));
	COVER_instance->rmatval=(double*) calloc(COVER_instance->nzcnt,sizeof(double));

	for(i=0; i<COVER_instance->cover_size; i++)
	{
		COVER_instance->rmatval[i]=weights[i];
		COVER_instance->rmatind[i]=i;
	}

	COVER_instance->rmatbeg[0]=0;

	COVER_instance->status=CPXaddrows(COVER_instance->env_sep[index],COVER_instance->lp_sep[index],0,COVER_instance->rcnt,COVER_instance->nzcnt,COVER_instance->rhs,COVER_instance->sense,COVER_instance->rmatbeg,COVER_instance->rmatind,COVER_instance->rmatval,NULL,NULL);
	if(COVER_instance->status!=0)
	{
		printf("error in CPXaddrows\n");
		exit(-1);
	}

	free(COVER_instance->rmatbeg);
	free(COVER_instance->rmatval);
	free(COVER_instance->rmatind);
	free(COVER_instance->rhs);
	free(COVER_instance->sense);
	/////////////////////////////////////////////////////////////////////////////////////////////////////////


	// * Set  tolerance *
	COVER_instance->status = CPXsetdblparam (COVER_instance->env_sep[index], CPX_PARAM_EPAGAP, 0.0);
	if (COVER_instance->status)
	{
		printf ("error for CPX_PARAM_EPAGAP\n");
	}

	// * Set  tolerance *
	COVER_instance->status = CPXsetdblparam (COVER_instance->env_sep[index], CPX_PARAM_EPGAP, 0.0);
	if (COVER_instance->status)
	{
		printf ("error for CPX_PARAM_EPGAP\n");
	}

	// * Set video output *
	COVER_instance->status = CPXsetintparam (COVER_instance->env_sep[index], CPX_PARAM_SCRIND, CPX_OFF);
	if (COVER_instance->status)
	{
		printf ("error for CPX_PARAM_SCRIND\n");
	}


	//	// * Set number of CPU*
	//
	COVER_instance->status = CPXsetintparam (COVER_instance->env_sep[index], CPX_PARAM_THREADS, 1);
	if (COVER_instance->status)
	{
		printf ("error for CPX_PARAM_THREADS\n");
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	// * writing the created ILP model on a file *
	//	COVER_instance->status=CPXwriteprob(COVER_instance->env_sep[index],COVER_instance->lp_sep[index],"sep.lp",NULL);
	//	if(COVER_instance->status!=0) {
	//		printf("error in CPXwriteprob\n");
	//		exit(-1);
	//	}
	/////////////////////////////////////////////////////////////////////////////////////////////////////////

}



/***********************************************************************************/
double solve_cplex_cover(data_COVER *COVER_instance,int index)
/***********************************************************************************/
{

	int i,j;


	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	COVER_instance->rmatind=(int*) calloc(COVER_instance->cover_size,sizeof(int));
	COVER_instance->rmatval=(double*) calloc(COVER_instance->cover_size,sizeof(double));
	for(i=0; i<COVER_instance->cover_size; i++)
	{
		COVER_instance->rmatind[i]=i;
		COVER_instance->rmatval[i]=-COVER_instance->profit_sep[i];
	}
	COVER_instance->status = CPXchgobj (COVER_instance->env_sep[index],COVER_instance->lp_sep[index], COVER_instance->cover_size, COVER_instance->rmatind, COVER_instance->rmatval);
	free(COVER_instance->rmatind);
	free(COVER_instance->rmatval);
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	///////////////////////////////////////////////////////////////////////////////////
	// * solving the MIP model
	clock_t time_start=clock();


	COVER_instance->status=CPXmipopt(COVER_instance->env_sep[index],COVER_instance->lp_sep[index]);
	if(COVER_instance->status!=0)
	{
		printf("error in CPXmipopt\n");
		exit(-1);
	}

	clock_t time_end=clock();
	double solution_time=(double)(time_end-time_start)/(double)CLOCKS_PER_SEC;

	//printf("\n\nsolution_time ->\t\%f",solution_time);

	///////////////////////////////////////////////////////////////////////////////////


	// * getting the solution

	COVER_instance->status=CPXgetmipx(COVER_instance->env_sep[index],COVER_instance->lp_sep[index],COVER_instance->cover,0,COVER_instance->cover_size-1);
	if(COVER_instance->status!=0)
	{
		printf("error in CPXgetmipx\n");
		exit(-1);
	}

	COVER_instance->status=CPXgetmipobjval(COVER_instance->env_sep[index],COVER_instance->lp_sep[index],&(COVER_instance->objval_cover));
	if(COVER_instance->status!=0)
	{
		printf("error in CPXgetmipobjval\n");
		exit(-1);
	}

	return COVER_instance->objval_cover+1;

}







/***********************************************************************************/
void free_cplex_cover(data_COVER *COVER_instance, int index)
/***********************************************************************************/
{


	COVER_instance->status=CPXfreeprob(COVER_instance->env_sep[index],&(COVER_instance->lp_sep[index]));
	if(COVER_instance->status!=0) {printf("error in CPXfreeprob\n");exit(-1);}

	COVER_instance->status=CPXcloseCPLEX(&(COVER_instance->env_sep[index]));
	if(COVER_instance->status!=0) {printf("error in CPXcloseCPLEX\n");exit(-1);}

}





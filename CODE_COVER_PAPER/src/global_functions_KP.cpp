

#include "global_functions_KP.h"

#define print_ist_features


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
void read_instance_file_KP(data *KP_instance)
/*********************************/
{


	int j;



	ifstream in(KP_instance->istname);
	if(!in)
	{
		ofstream err("Error.log",ios::app);
		cout << "File could not be opened. " << endl;
		err << "File could not be opened. " << endl;
		err.close();
		exit(1);
	}


	//cout << "INSTANCE->\t" << istname << endl;

	in >> KP_instance->item_number;
	in >> KP_instance->capacity;

	cout << "\n-->item_number\t"<<KP_instance->item_number<<endl;
	cout << "\n-->capacity\t"<<KP_instance->capacity<<endl;

	KP_instance->weights=(int*) calloc(KP_instance->item_number,sizeof(int));
	KP_instance->profits=(double*) calloc(KP_instance->item_number,sizeof(double));



	for(int i=0; i<KP_instance->item_number; i++)
	{
		in >>KP_instance->profits[i];
		in >>KP_instance->weights[i];
	}


	for(j=0; j<KP_instance->item_number; j++)
	{
		if(KP_instance->weights[j]>=KP_instance->capacity){KP_instance->weights[j]=KP_instance->capacity;}
	}


#ifdef print_ist_features

	cout << "Weights\t\t";
	for(j=0; j<KP_instance->item_number; j++)
	{
		cout << KP_instance->weights[j] << "\t";
	}
	cout << "\n";
	cout << "Profits\t\t";
	for(j=0; j<KP_instance->item_number; j++)
	{
		cout << KP_instance->profits[j]<< "\t";
	}
	cout << "\n";


#endif

	in.close();

	cout << "Instance read!\n";
	///////////////////////////////////////////////////////////

}


/*********************************/
void generate_instance_KP(data *KP_instance,int number_of_items,int perc_cap,int category,int R)
/*********************************/
{


	int j;


	//cout << "INSTANCE->\t" << istname << endl;

	KP_instance->item_number=number_of_items;

	KP_instance->weights=(int*) calloc(KP_instance->item_number,sizeof(int));
	KP_instance->profits=(double*) calloc(KP_instance->item_number,sizeof(double));

	long sum_weight=0;

	if(category==1)
	{

		cout << "Uncorrelated\n";

		for(int i=0; i<KP_instance->item_number; i++)
		{

			KP_instance->profits[i]=randomBETWEEN(1,R);

			KP_instance->weights[i]=randomBETWEEN(1,R);

			sum_weight+=KP_instance->weights[i];
		}
	}



	if(category==2)
	{

		cout <<  "\nWeakly correlated\n";

		for(int i=0; i<KP_instance->item_number; i++)
		{

			int dummy=randomBETWEEN(1,R);

			KP_instance->weights[i]=dummy;


			int LB_P=max( 1.0 , KP_instance->weights[i] - (R / 10.0) );

			int UB_P= KP_instance->weights[i] + (R / 10.0);

			KP_instance->profits[i]=randomBETWEEN(LB_P,UB_P);

			sum_weight+=KP_instance->weights[i];
		}
	}


	if(category==3)
	{

		cout <<  "\nStrongly correlated\n";

		for(int i=0; i<KP_instance->item_number; i++)
		{

			int dummy=randomBETWEEN(1,R);

			KP_instance->weights[i]=dummy;

			KP_instance->profits[i]=dummy + (R / 10);

			sum_weight+=KP_instance->weights[i];

		}
	}

	if(category==4)
	{
		cout <<  "\nInverse strongly correlated\n";

		for(int i=0; i<KP_instance->item_number; i++)
		{

			int dummy=randomBETWEEN(1,R);

			KP_instance->profits[i]=dummy;

			KP_instance->weights[i]= dummy + (R / 10);

			sum_weight+=KP_instance->weights[i];
		}

	}

	if(category==5)
	{
		cout <<  "\nAlmost strongly correlated\n";

		for(int i=0; i<KP_instance->item_number; i++)
		{

			int dummy=randomBETWEEN(1,R);

			KP_instance->weights[i]=dummy;

			int LB_P= KP_instance->weights[i] + (R / 10.0) - (R / 500.0);

			int UB_P= KP_instance->weights[i] + (R / 10.0) + (R / 500.0);

			KP_instance->profits[i]=randomBETWEEN(LB_P,UB_P);


			sum_weight+=KP_instance->weights[i];
		}
	}


	if(category==6)
	{

		cout << "Subset-sum\n";

		for(int i=0; i<KP_instance->item_number; i++)
		{

			int dummy=randomBETWEEN(1,R);

			KP_instance->profits[i]=dummy;

			KP_instance->weights[i]=dummy;

			sum_weight+=KP_instance->weights[i];
		}
	}

	if(sum_weight==0)
	{
		cout << "ERROR IN GENRATION\n";
		exit(-1);
	}



	KP_instance->capacity= (int) (sum_weight * perc_cap/100.0);


	for(j=0; j<KP_instance->item_number; j++)
	{
		if(KP_instance->weights[j]>=KP_instance->capacity){KP_instance->weights[j]=KP_instance->capacity;}
	}

	cout << "\n-->item_number\t"<<KP_instance->item_number<<endl;
	cout << "\n-->capacity\t"<<KP_instance->capacity<< "\tsum_weight\t" << sum_weight  << "\t%" << KP_instance->capacity/sum_weight <<endl;



#ifdef print_ist_features

	cout << "Weights\t\t";
	for(j=0; j<KP_instance->item_number; j++)
	{
		cout << KP_instance->weights[j] << "\t";
	}
	cout << "\n";
	cout << "Profits\t\t";
	for(j=0; j<KP_instance->item_number; j++)
	{
		cout << KP_instance->profits[j]<< "\t";
	}
	cout << "\n";


#endif


	cout << "Instance read!\n";
	///////////////////////////////////////////////////////////

}

/*********************************/
void free_data_prob(data *KP_instance)
/*********************************/
{

	free(KP_instance->weights);
	free(KP_instance->profits);
}


/***********************************************************************************/
void kp_load_cplex(data *KP_instance)
/***********************************************************************************/
{

	int i,j;


	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// * setting the CPLEX environment

	//opening the environment
	KP_instance->env=CPXopenCPLEX(&(KP_instance->status));
	if(KP_instance->status!=0)
	{
		printf("cannot open CPLEX environment\n");
		exit(-1);
	}

	//opening the pointer to the problem
	KP_instance->lp=CPXcreateprob(KP_instance->env,&(KP_instance->status),"KP");
	if(KP_instance->status!=0)
	{
		printf("cannot create problem\n");
		exit(-1);
	}
	/////////////////////////////////////////////////////////////////////////////////////////////////////////

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// * creating the variables *
	KP_instance->ccnt=KP_instance->item_number;
	KP_instance->obj=(double*) calloc(KP_instance->ccnt,sizeof(double));
	KP_instance->lb=(double*) calloc(KP_instance->ccnt,sizeof(double));
	KP_instance->ub=(double*) calloc(KP_instance->ccnt,sizeof(double));
	KP_instance->xctype=(char*) calloc(KP_instance->ccnt,sizeof(char));


	KP_instance->colname=(char**) calloc(KP_instance->ccnt,sizeof(char*));
	for(i=0;i<KP_instance->ccnt;i++){KP_instance->colname[i]=(char*) calloc(2000,sizeof(char));}

	for(i=0; i<KP_instance->ccnt; i++)
	{
		KP_instance->obj[i]=KP_instance->profits[i];
		KP_instance->lb[i]=0.0;
		//		KP_instance->ub[i]=CPX_INFBOUND;
		KP_instance->ub[i]=1.0;
		KP_instance->xctype[i]='B';
		sprintf(KP_instance->colname[i], "x%d",i+1);
	}

	KP_instance->status=CPXnewcols(KP_instance->env,KP_instance->lp,KP_instance->ccnt,KP_instance->obj,KP_instance->lb,KP_instance->ub,KP_instance->xctype,KP_instance->colname);
	if(KP_instance->status!=0)
	{
		printf("error in CPXnewcols\n");
		exit(-1);
	}

	free(KP_instance->obj);
	free(KP_instance->lb);
	free(KP_instance->ub);
	free(KP_instance->xctype);

	for(i=0;i<KP_instance->ccnt;i++){free(KP_instance->colname[i]);}
	free(KP_instance->colname);


	// * setting the objective function in the maximization form
	CPXchgobjsen(KP_instance->env,KP_instance->lp,CPX_MAX);
	/////////////////////////////////////////////////////////////////////////////////////////////////////////

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// * creating the knapsack constraint *
	KP_instance->rcnt=1;
	KP_instance->nzcnt=KP_instance->item_number;
	KP_instance->rhs=(double*) calloc(KP_instance->rcnt,sizeof(double));
	KP_instance->sense=(char*) calloc(KP_instance->rcnt,sizeof(double));

	KP_instance->rhs[0]=KP_instance->capacity;
	KP_instance->sense[0]='L';


	KP_instance->rmatbeg=(int*) calloc(KP_instance->rcnt,sizeof(int));
	KP_instance->rmatind=(int*) calloc(KP_instance->nzcnt,sizeof(int));
	KP_instance->rmatval=(double*) calloc(KP_instance->nzcnt,sizeof(double));

	for(i=0; i<KP_instance->item_number; i++)
	{
		KP_instance->rmatval[i]=KP_instance->weights[i];
		KP_instance->rmatind[i]=i;
	}

	KP_instance->rmatbeg[0]=0;

	KP_instance->status=CPXaddrows(KP_instance->env,KP_instance->lp,0,KP_instance->rcnt,KP_instance->nzcnt,KP_instance->rhs,KP_instance->sense,KP_instance->rmatbeg,KP_instance->rmatind,KP_instance->rmatval,NULL,NULL);
	if(KP_instance->status!=0)
	{
		printf("error in CPXaddrows\n");
		exit(-1);
	}

	free(KP_instance->rmatbeg);
	free(KP_instance->rmatval);
	free(KP_instance->rmatind);
	free(KP_instance->rhs);
	free(KP_instance->sense);
	/////////////////////////////////////////////////////////////////////////////////////////////////////////


	//	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	// * writing the created ILP model on a file *
	//	KP_instance->status=CPXwriteprob(KP_instance->env,KP_instance->lp,"kp.lp",NULL);
	//	if(KP_instance->status!=0) {
	//		printf("error in CPXwriteprob\n");
	//		exit(-1);
	//	}
	//	/////////////////////////////////////////////////////////////////////////////////////////////////////////




}

/***********************************************************************************/
void kp_free_cplex(data *KP_instance)
/***********************************************************************************/
{


	KP_instance->status=CPXfreeprob(KP_instance->env,&(KP_instance->lp));
	if(KP_instance->status!=0) {printf("error in CPXfreeprob\n");exit(-1);}

	KP_instance->status=CPXcloseCPLEX(&(KP_instance->env));
	if(KP_instance->status!=0) {printf("error in CPXcloseCPLEX\n");exit(-1);}

}

/***********************************************************************************/
double kp_solve_cplex(data *KP_instance)
/***********************************************************************************/
{

	int i,j;

	// * Set video output *
	KP_instance->status = CPXsetintparam (KP_instance->env, CPX_PARAM_SCRIND, CPX_OFF);
	if (KP_instance->status)
	{
		printf ("error for CPX_PARAM_SCRIND\n");
	}


	//	// * Set  tolerance *
	//	KP_instance->status = CPXsetdblparam (KP_instance->env, CPX_PARAM_EPAGAP, 0.0);
	//	if (KP_instance->status)
	//	{
	//		printf ("error for CPX_PARAM_EPAGAP\n");
	//	}
	//
	//	// * Set  tolerance *
	//	KP_instance->status = CPXsetdblparam (KP_instance->env, CPX_PARAM_EPGAP, 0.0);
	//	if (KP_instance->status)
	//	{
	//		printf ("error for CPX_PARAM_EPGAP\n");
	//	}
	//
	//
	//	// * Set mip tolerances integrality *
	//	KP_instance->status = CPXsetdblparam (KP_instance->env, CPX_PARAM_EPINT, 0.0);
	//	if (KP_instance->status)
	//	{
	//		printf ("error for CPX_PARAM_EPINTP\n");
	//	}
	//
	//	// * Set Feasibility tolerance *
	//	KP_instance->status = CPXsetdblparam (KP_instance->env, CPX_PARAM_EPRHS, 1e-9);
	//	if (KP_instance->status)
	//	{
	//		printf ("error for CPX_PARAM_EPRHS\n");
	//	}
	//
	//	// * Set number of CPU*
	//
	//	KP_instance->status = CPXsetintparam (KP_instance->env, CPX_PARAM_THREADS, KP_instance->number_of_CPU);
	//	if (KP_instance->status)
	//	{
	//		printf ("error for CPX_PARAM_THREADS\n");
	//	}
	//
	// * Set time limit *

	KP_instance->status = CPXsetdblparam (KP_instance->env, CPX_PARAM_TILIM,10);
	if (KP_instance->status)
	{
		printf ("error for CPX_PARAM_TILIM\n");
	}


	//	// * Set number of CPU*
	//
	KP_instance->status = CPXsetintparam (KP_instance->env, CPX_PARAM_THREADS, 1);
	if (KP_instance->status)
	{
		printf ("error for CPX_PARAM_THREADS\n");
	}

	///////////////////////////////////////////////////////////////////////////////////
	// * solving the MIP model
	clock_t time_start=clock();


	KP_instance->status=CPXmipopt(KP_instance->env,KP_instance->lp);
	if(KP_instance->status!=0)
	{
		printf("error in CPXmipopt\n");
		exit(-1);
	}

	clock_t time_end=clock();
	double solution_time=(double)(time_end-time_start)/(double)CLOCKS_PER_SEC;

	printf("\n\nsolution_time ->\t\%f",solution_time);


	///////////////////////////////////////////////////////////////////////////////////


	// * getting the solution

	KP_instance->xx=(double*) calloc(KP_instance->item_number,sizeof(double));


	KP_instance->status=CPXgetmipx(KP_instance->env,KP_instance->lp,KP_instance->xx,0,KP_instance->item_number-1);
	if(KP_instance->status!=0)
	{
		printf("error in CPXgetmipx\n");
		exit(-1);
	}

	KP_instance->status=CPXgetmipobjval(KP_instance->env,KP_instance->lp,&(KP_instance->objval));
	if(KP_instance->status!=0)
	{
		printf("error in CPXgetmipobjval\n");
		exit(-1);
	}

	printf("\n\nMIP solution value ->\t\%f",KP_instance->objval);


	//	printf("\n\nSolution\n");
	//	for (i = 0; i < KP_instance->item_number; i++)
	//	{
	//		printf("item %d -> %d\n",i+1 ,(int)(KP_instance->x[i]+0.5));
	//	}
	//	printf("\n\n");

	KP_instance->status=CPXgetbestobjval(KP_instance->env,KP_instance->lp,&(KP_instance->bestobjval));
	if(KP_instance->status!=0)
	{
		printf("error in CPXgetbestobjval\n");
		exit(-1);
	}

	KP_instance->lpstat=CPXgetstat(KP_instance->env,KP_instance->lp);
	KP_instance->nodecount = CPXgetnodecnt(KP_instance->env, KP_instance->lp);


	///////////////////////////////////////////////////////////////////////////////////
	/* linear programming relaxation*/

	double solution_time_lp=0;

	KP_instance->status = CPXchgprobtype (KP_instance->env, KP_instance->lp, CPXPROB_LP);

	clock_t time_start_lp=clock();
	KP_instance->status=CPXlpopt(KP_instance->env,KP_instance->lp);
	if(KP_instance->status!=0)
	{
		printf("err_FILEor in CPXlpopt slave solve\n");
		exit(-1);
	}

	clock_t time_end_lp=clock();

	solution_time_lp=(double)(time_end_lp-time_start_lp)/(double)CLOCKS_PER_SEC;

	double cplex_lp;

	KP_instance->status=CPXgetobjval(KP_instance->env,KP_instance->lp,&cplex_lp);
	if(KP_instance->status!=0)
	{
		printf("error in CPXgetmipobjval\n");
		exit(-1);
	}

	printf("\n\nLP solution value ->\t\%f",cplex_lp);

	KP_instance->status=CPXgetx(KP_instance->env,KP_instance->lp,KP_instance->xx,0,KP_instance->item_number-1);
	if(KP_instance->status!=0)
	{
		printf("error in CPXgetmipx\n");
		exit(-1);
	}

	//	printf("\n\nSolution\n");
	//	for (i = 0; i < KP_instance->item_number; i++)
	//	{
	//		printf("item %d -> %.3f (p/w %.3f)\n",i+1,KP_instance->xx[i],KP_instance->profits[i]/KP_instance->weights[i]);
	//	}
	//	printf("\n");

	///////////////////////////////////////////////////////////////////////////////////

	KP_instance->cur_numcols=CPXgetnumcols(KP_instance->env,KP_instance->lp);
	KP_instance->cur_numrows=CPXgetnumrows(KP_instance->env,KP_instance->lp);

	//printf("\nnumcols\t%d\n",KP_instance->cur_numcols);
	//printf("\nnumrows\t%d\n",KP_instance->cur_numrows);

	FILE *out;
	out=fopen("info_kp.txt","a+");
	fprintf(out,"%s\t%f\t%f\t%f\t%f\t%d\t%d\t%f\t%f\n",KP_instance->istname,KP_instance->capacity,KP_instance->objval,KP_instance->bestobjval,solution_time,KP_instance->lpstat,KP_instance->nodecount,cplex_lp,solution_time_lp);
	fclose(out);

	free(KP_instance->xx);

	return KP_instance->objval;

}


/***********************************************************************************/
double kp_solve_cplex_LP(data *KP_instance,double *point)
/***********************************************************************************/
{

	int i,j;

	// * Set video output *
	KP_instance->status = CPXsetintparam (KP_instance->env, CPX_PARAM_SCRIND, CPX_OFF);
	if (KP_instance->status)
	{
		printf ("error for CPX_PARAM_SCRIND\n");
	}

	// * Set number of CPU*

	KP_instance->status = CPXsetintparam (KP_instance->env, CPX_PARAM_THREADS, 1);
	if (KP_instance->status)
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

	KP_instance->status = CPXchgprobtype (KP_instance->env, KP_instance->lp, CPXPROB_LP);

	clock_t time_start_lp=clock();
	KP_instance->status=CPXlpopt(KP_instance->env,KP_instance->lp);
	if(KP_instance->status!=0)
	{
		printf("err_FILEor in CPXlpopt slave solve\n");
		exit(-1);
	}

	clock_t time_end_lp=clock();

	solution_time_lp=(double)(time_end_lp-time_start_lp)/(double)CLOCKS_PER_SEC;

	double cplex_lp;

	KP_instance->status=CPXgetobjval(KP_instance->env,KP_instance->lp,&cplex_lp);
	if(KP_instance->status!=0)
	{
		printf("error in CPXgetmipobjval\n");
		exit(-1);
	}

	printf("\n\nLP solution value ->\t\%f",cplex_lp);

	KP_instance->status=CPXgetx(KP_instance->env,KP_instance->lp,point,0,KP_instance->item_number-1);
	if(KP_instance->status!=0)
	{
		printf("error in CPXgetmipx\n");
		exit(-1);
	}

	//		printf("\n\nSolution\n");
	//		for (i = 0; i < KP_instance->item_number; i++)
	//		{
	//			printf("item %d -> %.3f (p/w %.3f)\n",i+1,point[i],KP_instance->profits[i]/KP_instance->weights[i]);
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


/***********************************************************************************/
void add_cover(data *KP_instance,double *cover)
/***********************************************************************************/
{

	int i,j;

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// * creating the knapsack constraint *
	KP_instance->rcnt=1;
	KP_instance->nzcnt=KP_instance->item_number;
	KP_instance->rhs=(double*) calloc(KP_instance->rcnt,sizeof(double));
	KP_instance->sense=(char*) calloc(KP_instance->rcnt,sizeof(double));

	KP_instance->sense[0]='L';


	KP_instance->rmatbeg=(int*) calloc(KP_instance->rcnt,sizeof(int));
	KP_instance->rmatind=(int*) calloc(KP_instance->nzcnt,sizeof(int));
	KP_instance->rmatval=(double*) calloc(KP_instance->nzcnt,sizeof(double));

	int counter=0;
	for(i=0; i<KP_instance->item_number; i++)
	{
		if(cover[i]>0.5)
		{
			KP_instance->rmatval[i]=1.0;
			counter++;
		}
		else{
			KP_instance->rmatval[i]=0.0;
		}
		KP_instance->rmatind[i]=i;
	}

	KP_instance->rhs[0]=counter-1;


	KP_instance->rmatbeg[0]=0;

	KP_instance->status=CPXaddrows(KP_instance->env,KP_instance->lp,0,KP_instance->rcnt,KP_instance->nzcnt,KP_instance->rhs,KP_instance->sense,KP_instance->rmatbeg,KP_instance->rmatind,KP_instance->rmatval,NULL,NULL);
	if(KP_instance->status!=0)
	{
		printf("error in CPXaddrows\n");
		exit(-1);
	}

	free(KP_instance->rmatbeg);
	free(KP_instance->rmatval);
	free(KP_instance->rmatind);
	free(KP_instance->rhs);
	free(KP_instance->sense);
	/////////////////////////////////////////////////////////////////////////////////////////////////////////


	//	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	// * writing the created ILP model on a file *
	//	KP_instance->status=CPXwriteprob(KP_instance->env,KP_instance->lp,"sep_cover.lp",NULL);
	//	if(KP_instance->status!=0) {
	//		printf("error in CPXwriteprob\n");
	//		exit(-1);
	//	}
	//	/////////////////////////////////////////////////////////////////////////////////////////////////////////

}




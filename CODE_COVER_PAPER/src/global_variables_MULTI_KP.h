#ifndef VARIABLE_local_HEADER
#define VARIABLE_local_HEADER


#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <set>
#include <string>
#include <sstream>
#include <errno.h>
#include <cstring>
#include <cstdlib>

using namespace std;

#include </home/fabio/ILOG/CPLEX_Studio_AcademicResearch129/cplex/include/ilcplex/cplex.h>


typedef struct data_MULTI_KP
{


	////////////////////////////////INPUT//////////////////////////
	int number_of_items;
	int perc_cap;
	int category;
	int R;
	int seed;
	int algo;
	int elle;
	int order_items_DP;
	double conflict_density;
	int number_of_KPs;
	int max_cuts;
	////////////////////////////////////////////////////////////////

	char *instname;
	double timeLimit;


	int status,ccnt,rcnt,nzcnt,lpstat,nodecount,cur_numrows,cur_numcols,stat;
	int* rmatbeg,*rmatind,*cmatbeg, *cmatind;
	double *rmatval,*cmatval,*rngval,*x,*pi,*obj, *lb, *ub,*rhs,coef_p,objval,bestobjval;
	char *xctype,*sense;
	char **rowname,**colname;
	double solution_time;

	double cplex_lp;
	double solution_time_lp;

	CPXENVptr env_MULTI_KP;
	CPXLPptr  lp_MULTI_KP;

	//////////////////////////////////////////////
	int N_ITEMS;
	int N_KNAPSACKS;
	int **MATRIX_WEIGHTS;
	double *MULTI_PROFITS;
	int *RHS_CAPACITIES;
	//////////////////////////////////////////////

	///////////////////////////////////////////////////////////////////////////////
	//conflicts of the MULTI KP instance
	int n_conflicts;
	int **CONF_MATRIX;//1 if there is a conflict and 0 otherwise
	///////////////////////////////////////////////////////////////////////////////

} data_MULTI_KP;

#endif

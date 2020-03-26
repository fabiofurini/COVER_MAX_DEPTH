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
//#include </Users/fabiofurini/Applications/IBM/ILOG/CPLEX_Studio127/cplex/include/ilcplex/cplex.h>



typedef struct data
{
	char *istname;
	int algorithm;
	int number_of_CPU;
	int timeLimit;

	int item_number;
	double capacity;
	int *weights;
	double *profits;

	int status,ccnt,rcnt,nzcnt,lpstat,nodecount,cur_numrows,cur_numcols;
	int* rmatbeg,*rmatind,*cmatbeg, *cmatind;
	double *rmatval,*cmatval,*rngval,*xx,*pi,*obj, *lb, *ub,*rhs,coef_p,objval,bestobjval;
	char *xctype,*sense;
	char **rowname,**colname;

	CPXENVptr env;
	CPXLPptr  lp;

//	CPXENVptr env_sep,env_sep_bis;
//	CPXLPptr  lp_sep,lp_sep_bis;
//	double *profit_sep;
//	//double *cover;
//	double objval_cover;
//	int cap_sep;
//	//int *x_int;



	//double *sol_DP;
	//vector < int > order_DP;


//	long int * magic_sum;
//	double *profit_sep_magic;

//	int number_of_covers;

} data;

#endif

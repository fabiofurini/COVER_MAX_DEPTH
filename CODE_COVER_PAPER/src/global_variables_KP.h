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

	///////////////////////////////////////////////////////////////////////////////
	//conflicts of the KP instance
	int n_conflicts;
	int **CONF_MATRIX;//1 if there is a conflict and 0 otherwise
	///////////////////////////////////////////////////////////////////////////////


} data;

#endif

#ifndef VARIABLE_COVER_local_HEADER
#define VARIABLE_COVER_local_HEADER


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

typedef struct data_COVER
{

	double magic_constant;

	//////////////////////////////////////////////////////////////////////////////////////////
	CPXENVptr env_sep;
	CPXLPptr  lp_sep;

	int status,ccnt,rcnt,nzcnt,lpstat,nodecount,cur_numrows,cur_numcols;
	int* rmatbeg,*rmatind,*cmatbeg, *cmatind;
	double *rmatval,*cmatval,*rngval,*xx,*pi,*obj, *lb, *ub,*rhs,coef_p,objval,bestobjval;
	char *xctype,*sense;
	char **rowname,**colname;
	//////////////////////////////////////////////////////////////////////////////////////////

	double * point;

	int cover_size;


	double *profit_sep;
	int *weight_sep;

	int capacity;

	double *cover;
	double objval_cover;
	int cap_sep;


	double *sol_DP;
	vector < int > order_DP;
	int max_i;


	long int * magic_sum;
	double *profit_sep_magic;

	int number_of_covers;

} data_COVER;

#endif

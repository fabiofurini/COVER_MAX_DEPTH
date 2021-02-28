#ifndef glo_COVER_function_HEADER
#define glo_COVER_function_HEADER


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
#include <algorithm>    // std::sort
#include <vector>       // std::vector

#include "global_variables_MULTI_KP.h"

#include "DP.h"

#include "global_variables_COVERs.h"

struct valuesSTR
{
	int item;
	double score;
};

/***********************************************************************************/
int extend_cover(int cover_size,int cap,double *cover,int *weight);
/***********************************************************************************/

/***********************************************************************************/
double solve_DP_lex(data_COVER *COVER_instance,int index);
/***********************************************************************************/

/***********************************************************************************/
int make_maximal(data_COVER *COVER_instance,int cover_size,int cap_sep,double *cover,int *weight_sep,double *profit_sep,int index);
/***********************************************************************************/

/***********************************************************************************/
void load_cplex_cover(data_COVER *COVER_instance,int index,int *weights, int capacity);
/***********************************************************************************/

/***********************************************************************************/
double solve_cplex_cover(data_COVER *COVER_instance,int index);
/***********************************************************************************/

/***********************************************************************************/
void free_cplex_cover(data_COVER *COVER_instance, int index);
/***********************************************************************************/

/***********************************************************************************/
void  COVER_DATA_free(data_COVER *COVER_instance,int n_KP);
/***********************************************************************************/

/***********************************************************************************/
void  COVER_DATA_init(data_COVER *COVER_instance,int cover_size,double *profits,int **weights, int *capacity,int n_KP,int shuffle);
/***********************************************************************************/


#endif

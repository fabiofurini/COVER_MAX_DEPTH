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


#include "global_variables_KP.h"
#include "global_variables_COVER.h"
#include "DP.h"

struct valuesSTR
{
	int item;
	double score;
};

/***********************************************************************************/
int extend_cover(int cover_size,int cap,double *cover,int *weight);
/***********************************************************************************/

/***********************************************************************************/
int make_maximal(int cover_size,int cap_sep,double *cover,int *weight_sep,double *profit_sep);
/***********************************************************************************/

/***********************************************************************************/
void load_cplex_cover(data_COVER *COVER_instance);
/***********************************************************************************/

/***********************************************************************************/
double solve_cplex_cover(data_COVER *COVER_instance);
/***********************************************************************************/

/***********************************************************************************/
void free_cplex_cover(data_COVER *COVER_instance);
/***********************************************************************************/

/***********************************************************************************/
void  COVER_DATA_free(data_COVER *COVER_instance);
/***********************************************************************************/

/***********************************************************************************/
void  COVER_DATA_init(data_COVER *COVER_instance,int cover_size,double *profits,int *weights, int capacity,int shuffle);
/***********************************************************************************/

/***********************************************************************************/
double solve_DP_magic(data_COVER *COVER_instance);
/***********************************************************************************/

/***********************************************************************************/
double solve_DP_lex(data_COVER *COVER_instance);
/***********************************************************************************/



#endif

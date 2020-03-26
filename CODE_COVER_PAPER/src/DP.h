#ifndef DP_HEADER
#define DP_HEADER


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



/*****************************************************************************/
double DP_kp_RATIO_POWER_FOR_COVER(int n, int C, double* p, int * w,double *sol,int max_i,double elle,double constant);
/*****************************************************************************/

/*****************************************************************************/
void Load_SOLUTION(double *sol,int C, int n, int i_star,bool ***UU, int * w);
/*****************************************************************************/

/*****************************************************************************/
double DP_kp01_advanced(int n, int C, double* p, int * w,double *sol,bool MAXIMALITY,vector < int > order);
/*****************************************************************************/


#endif

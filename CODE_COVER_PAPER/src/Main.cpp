#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <float.h>
#include <vector>
#include <algorithm>
#include <set>

using namespace std;

#include "global_variables_KP.h"

#include "global_functions_KP.h"

#include "global_variables_COVER.h"

#include "global_functions_COVER.h"

#include "DP.h"

#define TOLL_SEP 0.000001

#define check_debug


/******************************************/
int main(int argc, char** argv)
/******************************************/
{


	int algo;
	int elle;

	data KP_instance;

	KP_instance.istname=(char*) calloc(2000,sizeof(char));

	int number_of_items;
	int perc_cap;
	int category;
	int R;
	int seed;

	//istname=new char[2000];

	cout << "argc\t" << argc <<endl;

	if (argc == 8)
	{

		number_of_items=atoi(argv[1]);
		perc_cap=atoi(argv[2]);
		category=atoi(argv[3]);
		R=atoi(argv[4]);
		seed=atoi(argv[5]);
		algo=atoi(argv[6]);
		elle=atoi(argv[7]);

		srand(seed);

	}
	else
	{
		cout << "ERROR NUMBER STANDARD PARAMETERS" << endl;
		exit(2);
	}

	cout << "-> number_of_items\t\t\t"<< number_of_items <<"\n";
	cout << "-> perc_cap\t\t"<< perc_cap <<"\n";
	cout << "-> category\t\t\t"<< category <<"\n";
	cout << "-> R\t\t"<< R <<"\n";
	cout << "-> seed\t\t\t"<< seed <<"\n";
	cout << "-> algo\t\t"<< algo <<"\n";
	cout << "-> elle\t\t"<< elle <<"\n";


	//strcpy(KP_instance.istname,"toy_catan.txt");
	//read_instance_file_KP(&KP_instance);

	strcpy(KP_instance.istname,"GENERATED");
	generate_instance_KP(&KP_instance,number_of_items,perc_cap,category,R);

	KP_instance.timeLimit=999;
	KP_instance.number_of_CPU=1;


	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	cout << "\n\n********************SOLVE THE KP PROBLEM\n";

	kp_load_cplex(&KP_instance);

	clock_t time_start_KP=clock();
	double val_KP_MIP=kp_solve_cplex(&KP_instance);
	clock_t time_end_KP=clock();
	double solution_time_KP=(double)(time_end_KP-time_start_KP)/(double)CLOCKS_PER_SEC;

	kp_free_cplex(&KP_instance);
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	cout << "\n\n********************OPTIMIZE OVER THE COVER CLOSURE\n";


	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	data_COVER COVER_instance;

	COVER_DATA_init(&COVER_instance,KP_instance.item_number,KP_instance.weights,KP_instance.capacity);

	load_cplex_cover(&COVER_instance);

	COVER_instance.number_of_covers=0;
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	cout << "max_i\t" << COVER_instance.max_i<< endl;
	cout << "cap_sep\t" << COVER_instance.cap_sep<< endl;


	kp_load_cplex(&KP_instance);


	double LP_CURRENT=-1;
	double LP_INIT=-1;

	clock_t time_start=clock();

	while(1)
	{

		LP_CURRENT=kp_solve_cplex_LP(&KP_instance,COVER_instance.point);

		if(COVER_instance.number_of_covers==0)
		{
			LP_INIT=LP_CURRENT;
		}

		COVER_instance.magic_constant=0;

		int item_profit_frac=0;
		int item_profit_zero=0;
		int item_profit_uno=0;

		for(int i=0; i<KP_instance.item_number; i++)
		{
			COVER_instance.profit_sep[i]=1-COVER_instance.point[i];
			COVER_instance.magic_constant+=COVER_instance.point[i];

			if(COVER_instance.profit_sep[i]>TOLL_SEP && COVER_instance.profit_sep[i]< 1 -TOLL_SEP)
			{
				item_profit_frac++;
			}
			if(COVER_instance.profit_sep[i]<TOLL_SEP)
			{
				item_profit_zero++;
			}
			if(COVER_instance.profit_sep[i]> 1 -TOLL_SEP)
			{
				item_profit_uno++;
			}
		}

		cout << "\nitem_profit_frac\t" << item_profit_frac << "\t total item \t" << KP_instance.item_number << endl;
		cout << "item_profit_zero\t" << item_profit_zero << "\t total item \t" << KP_instance.item_number << endl;
		cout << "item_profit_uno\t" << item_profit_uno << "\t total item \t" << KP_instance.item_number << endl;
		cout << item_profit_frac+item_profit_zero+item_profit_uno << endl;


		double cover_val=-1;


		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		if(algo==1)
		{
			cover_val=solve_cplex_cover(&COVER_instance);


			if(cover_val > TOLL_SEP)
			{

				add_cover(&KP_instance,COVER_instance.cover);

				COVER_instance.number_of_covers++;
			}
			else
			{
				break;
			}
		}
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////


		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		if(algo==3)
		{

			cout << "\n\nCONSTANT\t" <<(1-KP_instance.item_number+COVER_instance.magic_constant) << endl <<endl;

			//ATTENTION THAT IT LOAD THE COMPLEMENT OF THE COVER!
			cover_val=DP_kp_RATIO_POWER_FOR_COVER(KP_instance.item_number,COVER_instance.cap_sep,COVER_instance.profit_sep,KP_instance.weights,COVER_instance.cover,COVER_instance.max_i,(double)elle,(1-KP_instance.item_number+COVER_instance.magic_constant));


			for (int i = 0; i < KP_instance.item_number; i++)
			{

				//printf("q=%.5f\tx=%.0f\t\tw=%d\ty=%.5f\n", 1-COVER_instance.point[i],COVER_instance.cover[i],COVER_instance.weight_sep[i],COVER_instance.point[i]);

				COVER_instance.cover[i]=1-COVER_instance.cover[i];
			}

			double val_check=cover_val;

			if(cover_val>TOLL_SEP)
			{

				add_cover(&KP_instance,COVER_instance.cover);

				COVER_instance.number_of_covers++;
			}
			else
			{

#ifdef check_debug
				cover_val=solve_cplex_cover(&COVER_instance);

				if(cover_val > TOLL_SEP)
				{

					cout << "\n\nHORROR\t" << cover_val << "\tDP--->\t" << val_check<< "\t algo \t"<< algo<< endl;

					cin.get();

					exit(-1);
				}
#endif

				break;
			}


		}
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////


		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		if(algo==4)
		{

			cover_val=DP_kp01_advanced(COVER_instance.cover_size,COVER_instance.cap_sep, COVER_instance.profit_sep,COVER_instance.weight_sep,COVER_instance.cover,false,COVER_instance.order_DP);


			cover_val= cover_val + COVER_instance.magic_constant - COVER_instance.cover_size +1;

			for (int i = 0; i < COVER_instance.cover_size; i++)
			{
				COVER_instance.cover[i]=1-COVER_instance.cover[i];
			}


			double val_check=cover_val;

			if(cover_val > TOLL_SEP)
			{

				add_cover(&KP_instance,COVER_instance.cover);

				COVER_instance.number_of_covers++;
			}
			else
			{

#ifdef check_debug
				cover_val=solve_cplex_cover(&COVER_instance);

				if(cover_val > TOLL_SEP)
				{

					cout << "\n\nHORROR\t" << cover_val << "\tDP--->\t" << val_check<< "\t algo \t"<< algo<< endl;

					cin.get();

					exit(-1);
				}
#endif

				break;
			}


		}
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////


		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		if(algo==14)
		{

			cover_val=DP_kp01_advanced(COVER_instance.cover_size,COVER_instance.cap_sep, COVER_instance.profit_sep,COVER_instance.weight_sep,COVER_instance.cover,false,COVER_instance.order_DP);

			make_maximal(COVER_instance.cover_size,COVER_instance.cap_sep,COVER_instance.cover,COVER_instance.weight_sep,COVER_instance.profit_sep);

			cover_val= cover_val + COVER_instance.magic_constant - COVER_instance.cover_size +1;

			for (int i = 0; i < COVER_instance.cover_size; i++)
			{
				COVER_instance.cover[i]=1-COVER_instance.cover[i];
			}


			double val_check=cover_val;

			if(cover_val > TOLL_SEP)
			{

				add_cover(&KP_instance,COVER_instance.cover);

				COVER_instance.number_of_covers++;
			}
			else
			{

#ifdef check_debug
				cover_val=solve_cplex_cover(&COVER_instance);

				if(cover_val > TOLL_SEP)
				{

					cout << "\n\nHORROR\t" << cover_val << "\tDP--->\t" << val_check<< "\t algo \t"<< algo<< endl;

					cin.get();

					exit(-1);
				}
#endif

				break;
			}


		}
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////



		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		if(algo==5)
		{

			cover_val=DP_kp01_advanced(COVER_instance.cover_size,COVER_instance.cap_sep, COVER_instance.profit_sep,COVER_instance.weight_sep,COVER_instance.cover,true,COVER_instance.order_DP);


			cover_val= cover_val + COVER_instance.magic_constant - COVER_instance.cover_size +1;

			for (int i = 0; i < COVER_instance.cover_size; i++)
			{
				COVER_instance.cover[i]=1-COVER_instance.cover[i];
			}


			double val_check=cover_val;

			if(cover_val > TOLL_SEP)
			{

				add_cover(&KP_instance,COVER_instance.cover);

				COVER_instance.number_of_covers++;
			}
			else
			{

#ifdef check_debug
				cover_val=solve_cplex_cover(&COVER_instance);

				if(cover_val > TOLL_SEP)
				{
					cout << "\n\nHORROR\t" << cover_val << "\tDP--->\t" << val_check<< "\t algo \t"<< algo<< endl;

					cin.get();

					exit(-1);
				}
#endif

				break;
			}


		}
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////


		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		if(algo==6)
		{

			cover_val=solve_DP_magic(&COVER_instance);

			double val_check=cover_val;

			if(cover_val<1-TOLL_SEP)
			{

				add_cover(&KP_instance,COVER_instance.cover);

				COVER_instance.number_of_covers++;
			}
			else
			{

#ifdef check_debug
				cover_val=solve_cplex_cover(&COVER_instance);

				if(cover_val > TOLL_SEP)
				{
					cout << "\n\nHORROR\t" << cover_val << "\tDP--->\t" << val_check<< endl;

					cin.get();

					exit(-1);
				}
#endif

				break;
			}


		}
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		if(algo==7)
		{

			cover_val=solve_DP_lex(&COVER_instance);

			double val_check=cover_val;

			if(cover_val<1-TOLL_SEP)
			{

				add_cover(&KP_instance,COVER_instance.cover);

				COVER_instance.number_of_covers++;
			}
			else
			{

#ifdef check_debug
				cover_val=solve_cplex_cover(&COVER_instance);

				if(cover_val > TOLL_SEP)
				{
					cout << "\n\nHORROR\t" << cover_val << "\tDP--->\t" << val_check<< endl;

					cin.get();

					exit(-1);
				}
#endif

				break;
			}


		}
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////



		if(cover_val==-1){break;}

	}

	clock_t time_end=clock();
	double solution_time=(double)(time_end-time_start)/(double)CLOCKS_PER_SEC;


	kp_free_cplex(&KP_instance);

	free_cplex_cover(&COVER_instance);

	cout << "\n\nGENERATED\t" << COVER_instance.number_of_covers << "\tCOVERS" <<endl;
	cout << "TIME\t" << solution_time <<endl;
	cout << fixed << "LP_CURRENT\t" << LP_CURRENT << endl;
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	//	/////////////////////////////////////////////////////////////////////////////
	ofstream compact_file;
	compact_file.open("info_Exensive.txt", ios::app);
	compact_file << fixed

			<< "GENERATED" << "\t"

			<< algo << "\t"
			<< elle << "\t"

			<< KP_instance.item_number << "\t"
			<< KP_instance.capacity << "\t"

			<< LP_CURRENT << "\t"
			<< COVER_instance.number_of_covers << "\t"
			<< solution_time << "\t"

			<<  val_KP_MIP << "\t"
			<<  solution_time_KP << "\t"

			<<   LP_INIT << "\t"

			<<    COVER_instance.max_i << "\t"
			<< COVER_instance.cap_sep << "\t"


			<<  number_of_items << "\t"
			<< 	perc_cap<< "\t"
			<<  category<< "\t"
			<< 	R << "\t"
			<<  seed<< "\t"

			<< endl;
	compact_file.close();
	//	/////////////////////////////////////////////////////////////////////////////


	free(KP_instance.istname);

	//////////////////
	free_data_prob(&KP_instance);
	//////////////////

	COVER_DATA_free(&COVER_instance);

	printf("\nDONE!\n");

	return 1;
}



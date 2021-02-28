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


#include "global_functions_MULTI_KP.h"
#include "global_variables_MULTI_KP.h"
#include "DP.h"
#include "global_functions_COVERs.h"


#define TOLL_SEP 0.000001

/******************************************/
int main(int argc, char** argv)
/******************************************/
{


	data_MULTI_KP MULTI_KP_instance;

	MULTI_KP_instance.instname=(char*) calloc(2000,sizeof(char));

	//	cout << "argc\t" << argc <<endl;
	if (argc == 12)
	{

		MULTI_KP_instance.number_of_items=atoi(argv[1]);
		MULTI_KP_instance.perc_cap=atoi(argv[2]);
		MULTI_KP_instance.category=atoi(argv[3]);
		MULTI_KP_instance.R=atoi(argv[4]);
		MULTI_KP_instance.seed=atoi(argv[5]);
		MULTI_KP_instance.algo=atoi(argv[6]);
		MULTI_KP_instance.elle=atoi(argv[7]);
		MULTI_KP_instance.order_items_DP=atoi(argv[8]);
		MULTI_KP_instance.conflict_density=atof(argv[9]);
		MULTI_KP_instance.number_of_KPs=atoi(argv[10]);
		MULTI_KP_instance.max_cuts=atoi(argv[11]);

		srand(MULTI_KP_instance.seed);

	}
	else
	{
		cout << "ERROR NUMBER STANDARD PARAMETERS" << endl;
		exit(2);
	}

	cout << "-> number_of_items\t\t\t"<<MULTI_KP_instance.number_of_items<<"\n";
	cout << "-> perc_cap\t\t\t"<<MULTI_KP_instance.perc_cap<<"\n";
	cout << "-> category\t\t\t"<<MULTI_KP_instance.category<<"\n";
	cout << "-> R\t\t\t"<<MULTI_KP_instance.R<<"\n";
	cout << "-> seed\t\t\t"<<MULTI_KP_instance.seed<<"\n";
	cout << "-> algo\t\t\t"<<MULTI_KP_instance.algo<<"\n";
	cout << "-> elle\t\t\t"<<MULTI_KP_instance.elle<<"\n";
	cout << "-> order_items_DP\t\t\t"<<MULTI_KP_instance.order_items_DP<<"\n";
	cout << "-> conflict_density\t\t\t"<<MULTI_KP_instance.conflict_density<<"\n";
	cout << "-> number_of_KPs\t\t\t"<<MULTI_KP_instance.number_of_KPs<<"\n";
	cout << "-> max_cuts\t\t\t"<<MULTI_KP_instance.max_cuts<<"\n";

	strcpy(MULTI_KP_instance.instname,"TEST");
	MULTI_KP_instance.timeLimit=10000;


	generate_instance_MULTI_KP(&MULTI_KP_instance);


	//	//TO READ AN INSTANCE
	//	MULTI_KP_instance.n_conflicts=read_instance_file_multi_kp(&MULTI_KP_instance,MULTI_KP_instance.conflict_density);


	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	cout << "\n\n********************SOLVE THE MULTI-KP PROBLEM\n";

	multi_kp_load_cplex(&MULTI_KP_instance);

	multi_kp_solve_cplex(&MULTI_KP_instance);

	multi_kp_free_cplex(&MULTI_KP_instance);

	//	exit(-1);
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	cout << "\n\n********************OPTIMIZE OVER THE COVER CLOSURE\n";


	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	data_COVER COVER_instance;

	COVER_DATA_init(&COVER_instance,MULTI_KP_instance.N_ITEMS,MULTI_KP_instance.MULTI_PROFITS,MULTI_KP_instance.MATRIX_WEIGHTS,MULTI_KP_instance.RHS_CAPACITIES,MULTI_KP_instance.N_KNAPSACKS,MULTI_KP_instance.order_items_DP);

	COVER_instance.env_sep=(CPXENVptr*) calloc(MULTI_KP_instance.N_KNAPSACKS,sizeof(CPXENVptr));
	COVER_instance.lp_sep=(CPXLPptr*) calloc(MULTI_KP_instance.N_KNAPSACKS,sizeof(CPXLPptr));

	for(int i=0; i<MULTI_KP_instance.N_KNAPSACKS; i++)
	{
		load_cplex_cover(&COVER_instance,i,MULTI_KP_instance.MATRIX_WEIGHTS[i],MULTI_KP_instance.RHS_CAPACITIES[i]);
	}

	COVER_instance.number_of_covers=0;
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	for(int i=0; i<MULTI_KP_instance.N_KNAPSACKS; i++)
	{
		cout << "max_i\t" << COVER_instance.max_i[i]<< endl;
		cout << "cap_sep\t" << COVER_instance.cap_sep[i]<< endl;
	}

	multi_kp_load_cplex(&MULTI_KP_instance);


	double LP_CURRENT=-1;
	double LP_INIT=-1;

	clock_t time_start=clock();

	double density_total=0;
	double density_min=1.0;
	double density_max=0.0;


	double viol_total=0;
	double viol_min=100000.0;
	double viol_max=0.0;

	double item_total=0;
	double item_min=MULTI_KP_instance.N_ITEMS;
	double item_max=0.0;

	double total_item_profit_frac=0;
	double total_item_profit_zero=0;
	double total_item_profit_uno=0;
	int total_iterations=0;

	double total_items_extended=0;

	while(1)
	{

		LP_CURRENT=multi_kp_solve_cplex_LP(&MULTI_KP_instance,COVER_instance.point);

		if(COVER_instance.number_of_covers==0)
		{
			LP_INIT=LP_CURRENT;
		}

		COVER_instance.magic_constant=0;
		for(int i=0; i<MULTI_KP_instance.N_ITEMS; i++)
		{
			COVER_instance.profit_sep[i]=1-COVER_instance.point[i];
			COVER_instance.magic_constant+=COVER_instance.point[i];
		}


		if(COVER_instance.number_of_covers>=MULTI_KP_instance.max_cuts)
		{
			break;
		}

		double cover_val=-1;

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		int item_profit_frac=0;
		int item_profit_zero=0;
		int item_profit_uno=0;

		for(int i=0; i<MULTI_KP_instance.N_ITEMS; i++)
		{

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

		//		cout << "\nitem_profit_frac\t" << item_profit_frac << "\t total item \t" << KP_instance.item_number << endl;
		//		cout << "item_profit_zero\t" << item_profit_zero << "\t total item \t" << KP_instance.item_number << endl;
		//		cout << "item_profit_uno\t" << item_profit_uno << "\t total item \t" << KP_instance.item_number << endl;
		//		cout << item_profit_frac+item_profit_zero+item_profit_uno << endl;

		total_item_profit_frac+=item_profit_frac;
		total_item_profit_zero+=item_profit_zero;
		total_item_profit_uno+=item_profit_uno;

		total_iterations++;
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////





		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		if(MULTI_KP_instance.algo==1)
		{
			bool NOT_ADDED=true;

			for(int i=0; i<MULTI_KP_instance.N_KNAPSACKS; i++)
			{
				cover_val=solve_cplex_cover(&COVER_instance,i);

				if(cover_val>TOLL_SEP)
				{


					////////////////////////////////////////////////////////////////////////////////////////////////////////////
					double counter=0;
					for (int i = 0; i < MULTI_KP_instance.N_ITEMS; i++){if(COVER_instance.cover[i]>0.5){counter++;}}
					double viol=0;
					for (int i = 0; i < MULTI_KP_instance.N_ITEMS; i++){if(COVER_instance.cover[i]>0.5){viol+=(COVER_instance.point[i]-1);}}viol++;
					cout << "cover_val\t" << cover_val << "\tnumber_of_covers\t"<< COVER_instance.number_of_covers << "\t density \t" << counter/MULTI_KP_instance.N_ITEMS << "\t item \t" << counter << "\t viol \t" << viol << endl;
					viol_total+=viol;
					if(viol_max<viol){viol_max=viol;}
					if(viol_min>viol){viol_min=viol;}
					density_total+=counter/MULTI_KP_instance.N_ITEMS;
					if(density_max<counter/MULTI_KP_instance.N_ITEMS){density_max=counter/MULTI_KP_instance.N_ITEMS;}
					if(density_min>counter/MULTI_KP_instance.N_ITEMS){density_min=counter/MULTI_KP_instance.N_ITEMS;}
					item_total+=counter;
					if(item_max<counter){item_max=counter;}
					if(item_min>counter){item_min=counter;}
					////////////////////////////////////////////////////////////////////////////////////////////////////////////


					add_cover(&MULTI_KP_instance,COVER_instance.cover);

					NOT_ADDED=false;

					COVER_instance.number_of_covers++;
				}

			}

			if(NOT_ADDED)
			{
				break;
			}
		}
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////



		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		if(MULTI_KP_instance.algo==3)
		{

			bool NOT_ADDED=true;

			//cout << "\n\nCONSTANT\t" <<(1-MULTI_KP_instance.N_ITEMS+COVER_instance.magic_constant) << endl <<endl;


			for(int i=0; i<MULTI_KP_instance.N_KNAPSACKS; i++)
			{
				//ATTENTION THAT IT LOAD THE COMPLEMENT OF THE COVER!

				cover_val=DP_kp_RATIO_POWER_FOR_COVER(MULTI_KP_instance.N_ITEMS,COVER_instance.cap_sep[i],COVER_instance.profit_sep,MULTI_KP_instance.MATRIX_WEIGHTS[i],COVER_instance.cover,COVER_instance.max_i[i],(double)MULTI_KP_instance.elle,(1-MULTI_KP_instance.N_ITEMS+COVER_instance.magic_constant),COVER_instance.order_DP[i]);

				for (int i = 0; i < MULTI_KP_instance.N_ITEMS; i++)
				{
					COVER_instance.cover[i]=1-COVER_instance.cover[i];

				}

				if(cover_val>TOLL_SEP)
				{

					////////////////////////////////////////////////////////////////////////////////////////////////////////////
					double counter=0;
					for (int i = 0; i < MULTI_KP_instance.N_ITEMS; i++){if(COVER_instance.cover[i]>0.5){counter++;}}
					double viol=0;
					for (int i = 0; i < MULTI_KP_instance.N_ITEMS; i++){if(COVER_instance.cover[i]>0.5){viol+=(COVER_instance.point[i]-1);}}viol++;
					cout << "cover_val\t" << cover_val << "\tnumber_of_covers\t"<< COVER_instance.number_of_covers << "\t density \t" << counter/MULTI_KP_instance.N_ITEMS << "\t item \t" << counter << "\t viol \t" << viol << endl;
					viol_total+=viol;
					if(viol_max<viol){viol_max=viol;}
					if(viol_min>viol){viol_min=viol;}
					density_total+=counter/MULTI_KP_instance.N_ITEMS;
					if(density_max<counter/MULTI_KP_instance.N_ITEMS){density_max=counter/MULTI_KP_instance.N_ITEMS;}
					if(density_min>counter/MULTI_KP_instance.N_ITEMS){density_min=counter/MULTI_KP_instance.N_ITEMS;}
					item_total+=counter;
					if(item_max<counter){item_max=counter;}
					if(item_min>counter){item_min=counter;}
					////////////////////////////////////////////////////////////////////////////////////////////////////////////


					add_cover(&MULTI_KP_instance,COVER_instance.cover);

					COVER_instance.number_of_covers++;

					NOT_ADDED=false;
				}

			}

			if(NOT_ADDED)
			{
				break;
			}

		}
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////



		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		if(MULTI_KP_instance.algo==33)
		{

			bool NOT_ADDED=true;

			//cout << "\n\nCONSTANT\t" <<(1-MULTI_KP_instance.N_ITEMS+COVER_instance.magic_constant) << endl <<endl;


			for(int i=0; i<MULTI_KP_instance.N_KNAPSACKS; i++)
			{
				//ATTENTION THAT IT LOAD THE COMPLEMENT OF THE COVER!
				cover_val=DP_kp_RATIO_POWER_FOR_COVER(MULTI_KP_instance.N_ITEMS,COVER_instance.cap_sep[i],COVER_instance.profit_sep,MULTI_KP_instance.MATRIX_WEIGHTS[i],COVER_instance.cover,COVER_instance.max_i[i],(double)MULTI_KP_instance.elle,(1-MULTI_KP_instance.N_ITEMS+COVER_instance.magic_constant),COVER_instance.order_DP[i]);

				//cout << "\nval_DP\t" << cover_val << endl;

				for (int i = 0; i < MULTI_KP_instance.N_ITEMS; i++)
				{
					COVER_instance.cover[i]=1-COVER_instance.cover[i];

				}


				if(cover_val>TOLL_SEP)
				{

					////////////////////////////////////////////////////////////////////////////////////////////////////////////
					double counter=0;
					for (int i = 0; i < MULTI_KP_instance.N_ITEMS; i++){if(COVER_instance.cover[i]>0.5){counter++;}}
					double viol=0;
					for (int i = 0; i < MULTI_KP_instance.N_ITEMS; i++){if(COVER_instance.cover[i]>0.5){viol+=(COVER_instance.point[i]-1);}}viol++;
					cout << "cover_val (EXTENDED)\t" << cover_val << "\tnumber_of_covers\t"<< COVER_instance.number_of_covers << "\t density \t" << counter/MULTI_KP_instance.N_ITEMS<< "\t item \t" << counter << "\tviol\t " << viol<< endl;
					viol_total+=viol;
					if(viol_max<viol){viol_max=viol;}
					if(viol_min>viol){viol_min=viol;}
					density_total+=counter/MULTI_KP_instance.N_ITEMS;
					if(density_max<counter/MULTI_KP_instance.N_ITEMS){density_max=counter/MULTI_KP_instance.N_ITEMS;}
					if(density_min>counter/MULTI_KP_instance.N_ITEMS){density_min=counter/MULTI_KP_instance.N_ITEMS;}
					item_total+=counter;
					if(item_max<counter){item_max=counter;}
					if(item_min>counter){item_min=counter;}
					////////////////////////////////////////////////////////////////////////////////////////////////////////////


					/////////////////////////////////////////////////////////////////////////////////////////////////////////
					int eureka=extend_cover(MULTI_KP_instance.N_ITEMS,MULTI_KP_instance.RHS_CAPACITIES[i],COVER_instance.cover,MULTI_KP_instance.MATRIX_WEIGHTS[i]);
					cout << "eureka\t" << eureka << endl;
					/////////////////////////////////////////////////////////////////////////////////////////////////////////

					add_cover_extended(&MULTI_KP_instance,COVER_instance.cover,counter);

					total_items_extended+=eureka;

					COVER_instance.number_of_covers++;


					NOT_ADDED=false;
				}

			}

			if(NOT_ADDED)
			{
				break;
			}

		}
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////



		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		if(MULTI_KP_instance.algo==4)
		{

			bool NOT_ADDED=true;

			for(int i=0; i<MULTI_KP_instance.N_KNAPSACKS; i++)
			{

				cover_val=DP_kp01_advanced(COVER_instance.cover_size,COVER_instance.cap_sep[i], COVER_instance.profit_sep,COVER_instance.weight_sep[i],COVER_instance.cover,false,COVER_instance.order_DP[i]);

				cover_val= cover_val + COVER_instance.magic_constant - COVER_instance.cover_size +1;

				for (int i = 0; i < COVER_instance.cover_size; i++)
				{
					COVER_instance.cover[i]=1-COVER_instance.cover[i];
				}


				if(cover_val > TOLL_SEP)
				{


					////////////////////////////////////////////////////////////////////////////////////////////////////////////
					double counter=0;
					for (int i = 0; i < MULTI_KP_instance.N_ITEMS; i++){if(COVER_instance.cover[i]>0.5){counter++;}}
					double viol=0;
					for (int i = 0; i < MULTI_KP_instance.N_ITEMS; i++){if(COVER_instance.cover[i]>0.5){viol+=(COVER_instance.point[i]-1);}}viol++;
					cout << "cover_val\t" << cover_val << "\tnumber_of_covers\t"<< COVER_instance.number_of_covers << "\t density \t" << counter/MULTI_KP_instance.N_ITEMS << "\t item \t" << counter << "\t viol \t" << viol << endl;
					viol_total+=viol;
					if(viol_max<viol){viol_max=viol;}
					if(viol_min>viol){viol_min=viol;}
					density_total+=counter/MULTI_KP_instance.N_ITEMS;
					if(density_max<counter/MULTI_KP_instance.N_ITEMS){density_max=counter/MULTI_KP_instance.N_ITEMS;}
					if(density_min>counter/MULTI_KP_instance.N_ITEMS){density_min=counter/MULTI_KP_instance.N_ITEMS;}
					item_total+=counter;
					if(item_max<counter){item_max=counter;}
					if(item_min>counter){item_min=counter;}
					////////////////////////////////////////////////////////////////////////////////////////////////////////////


					add_cover(&MULTI_KP_instance,COVER_instance.cover);

					COVER_instance.number_of_covers++;

					NOT_ADDED=false;
				}

			}

			if(NOT_ADDED)
			{
				break;
			}


		}
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////



		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		if(MULTI_KP_instance.algo==44)
		{

			bool NOT_ADDED=true;

			for(int i=0; i<MULTI_KP_instance.N_KNAPSACKS; i++)
			{

				cover_val=DP_kp01_advanced(COVER_instance.cover_size,COVER_instance.cap_sep[i], COVER_instance.profit_sep,COVER_instance.weight_sep[i],COVER_instance.cover,false,COVER_instance.order_DP[i]);

				cover_val= cover_val + COVER_instance.magic_constant - COVER_instance.cover_size +1;

				for (int i = 0; i < COVER_instance.cover_size; i++)
				{
					COVER_instance.cover[i]=1-COVER_instance.cover[i];
				}


				if(cover_val > TOLL_SEP)
				{


					////////////////////////////////////////////////////////////////////////////////////////////////////////////
					double counter=0;
					for (int i = 0; i < MULTI_KP_instance.N_ITEMS; i++){if(COVER_instance.cover[i]>0.5){counter++;}}
					double viol=0;
					for (int i = 0; i < MULTI_KP_instance.N_ITEMS; i++){if(COVER_instance.cover[i]>0.5){viol+=(COVER_instance.point[i]-1);}}viol++;
					cout << "cover_val (EXTENDED)\t" << cover_val << "\tnumber_of_covers\t"<< COVER_instance.number_of_covers << "\t density \t" << counter/MULTI_KP_instance.N_ITEMS<< "\t item \t" << counter << "\tviol\t " << viol<< endl;
					viol_total+=viol;
					if(viol_max<viol){viol_max=viol;}
					if(viol_min>viol){viol_min=viol;}
					density_total+=counter/MULTI_KP_instance.N_ITEMS;
					if(density_max<counter/MULTI_KP_instance.N_ITEMS){density_max=counter/MULTI_KP_instance.N_ITEMS;}
					if(density_min>counter/MULTI_KP_instance.N_ITEMS){density_min=counter/MULTI_KP_instance.N_ITEMS;}
					item_total+=counter;
					if(item_max<counter){item_max=counter;}
					if(item_min>counter){item_min=counter;}
					////////////////////////////////////////////////////////////////////////////////////////////////////////////


					/////////////////////////////////////////////////////////////////////////////////////////////////////////
					int eureka=extend_cover(MULTI_KP_instance.N_ITEMS,MULTI_KP_instance.RHS_CAPACITIES[i],COVER_instance.cover,MULTI_KP_instance.MATRIX_WEIGHTS[i]);
					cout << "eureka\t" << eureka << endl;
					/////////////////////////////////////////////////////////////////////////////////////////////////////////

					add_cover_extended(&MULTI_KP_instance,COVER_instance.cover,counter);

					total_items_extended+=eureka;

					COVER_instance.number_of_covers++;

					NOT_ADDED=false;
				}

			}

			if(NOT_ADDED)
			{
				break;
			}


		}
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////



		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		if(MULTI_KP_instance.algo==5)
		{

			bool NOT_ADDED=true;

			for(int i=0; i<MULTI_KP_instance.N_KNAPSACKS; i++)
			{

				cover_val=DP_kp01_advanced(COVER_instance.cover_size,COVER_instance.cap_sep[i], COVER_instance.profit_sep,COVER_instance.weight_sep[i],COVER_instance.cover,true,COVER_instance.order_DP[i]);

				cover_val= cover_val + COVER_instance.magic_constant - COVER_instance.cover_size +1;

				for (int i = 0; i < COVER_instance.cover_size; i++)
				{
					COVER_instance.cover[i]=1-COVER_instance.cover[i];
				}


				if(cover_val > TOLL_SEP)
				{


					////////////////////////////////////////////////////////////////////////////////////////////////////////////
					double counter=0;
					for (int i = 0; i < MULTI_KP_instance.N_ITEMS; i++){if(COVER_instance.cover[i]>0.5){counter++;}}
					double viol=0;
					for (int i = 0; i < MULTI_KP_instance.N_ITEMS; i++){if(COVER_instance.cover[i]>0.5){viol+=(COVER_instance.point[i]-1);}}viol++;
					cout << "cover_val\t" << cover_val << "\tnumber_of_covers\t"<< COVER_instance.number_of_covers << "\t density \t" << counter/MULTI_KP_instance.N_ITEMS << "\t item \t" << counter << "\t viol \t" << viol << endl;
					viol_total+=viol;
					if(viol_max<viol){viol_max=viol;}
					if(viol_min>viol){viol_min=viol;}
					density_total+=counter/MULTI_KP_instance.N_ITEMS;
					if(density_max<counter/MULTI_KP_instance.N_ITEMS){density_max=counter/MULTI_KP_instance.N_ITEMS;}
					if(density_min>counter/MULTI_KP_instance.N_ITEMS){density_min=counter/MULTI_KP_instance.N_ITEMS;}
					item_total+=counter;
					if(item_max<counter){item_max=counter;}
					if(item_min>counter){item_min=counter;}
					////////////////////////////////////////////////////////////////////////////////////////////////////////////


					add_cover(&MULTI_KP_instance,COVER_instance.cover);

					COVER_instance.number_of_covers++;

					NOT_ADDED=false;
				}

			}

			if(NOT_ADDED)
			{
				break;
			}


		}
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////



		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		if(MULTI_KP_instance.algo==55)
		{

			bool NOT_ADDED=true;

			for(int i=0; i<MULTI_KP_instance.N_KNAPSACKS; i++)
			{

				cover_val=DP_kp01_advanced(COVER_instance.cover_size,COVER_instance.cap_sep[i], COVER_instance.profit_sep,COVER_instance.weight_sep[i],COVER_instance.cover,true,COVER_instance.order_DP[i]);

				cover_val= cover_val + COVER_instance.magic_constant - COVER_instance.cover_size +1;

				for (int i = 0; i < COVER_instance.cover_size; i++)
				{
					COVER_instance.cover[i]=1-COVER_instance.cover[i];
				}


				if(cover_val > TOLL_SEP)
				{


					////////////////////////////////////////////////////////////////////////////////////////////////////////////
					double counter=0;
					for (int i = 0; i < MULTI_KP_instance.N_ITEMS; i++){if(COVER_instance.cover[i]>0.5){counter++;}}
					double viol=0;
					for (int i = 0; i < MULTI_KP_instance.N_ITEMS; i++){if(COVER_instance.cover[i]>0.5){viol+=(COVER_instance.point[i]-1);}}viol++;
					cout << "cover_val (EXTENDED)\t" << cover_val << "\tnumber_of_covers\t"<< COVER_instance.number_of_covers << "\t density \t" << counter/MULTI_KP_instance.N_ITEMS<< "\t item \t" << counter << "\tviol\t " << viol<< endl;
					viol_total+=viol;
					if(viol_max<viol){viol_max=viol;}
					if(viol_min>viol){viol_min=viol;}
					density_total+=counter/MULTI_KP_instance.N_ITEMS;
					if(density_max<counter/MULTI_KP_instance.N_ITEMS){density_max=counter/MULTI_KP_instance.N_ITEMS;}
					if(density_min>counter/MULTI_KP_instance.N_ITEMS){density_min=counter/MULTI_KP_instance.N_ITEMS;}
					item_total+=counter;
					if(item_max<counter){item_max=counter;}
					if(item_min>counter){item_min=counter;}
					////////////////////////////////////////////////////////////////////////////////////////////////////////////


					/////////////////////////////////////////////////////////////////////////////////////////////////////////
					int eureka=extend_cover(MULTI_KP_instance.N_ITEMS,MULTI_KP_instance.RHS_CAPACITIES[i],COVER_instance.cover,MULTI_KP_instance.MATRIX_WEIGHTS[i]);
					cout << "eureka\t" << eureka << endl;
					/////////////////////////////////////////////////////////////////////////////////////////////////////////

					add_cover_extended(&MULTI_KP_instance,COVER_instance.cover,counter);

					total_items_extended+=eureka;

					COVER_instance.number_of_covers++;

					NOT_ADDED=false;
				}

			}

			if(NOT_ADDED)
			{
				break;
			}


		}
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////


		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		if(MULTI_KP_instance.algo==14)
		{

			bool NOT_ADDED=true;

			for(int i=0; i<MULTI_KP_instance.N_KNAPSACKS; i++)
			{

				cover_val=DP_kp01_advanced(COVER_instance.cover_size,COVER_instance.cap_sep[i], COVER_instance.profit_sep,COVER_instance.weight_sep[i],COVER_instance.cover,false,COVER_instance.order_DP[i]);

				make_maximal(&COVER_instance,COVER_instance.cover_size,COVER_instance.cap_sep[i],COVER_instance.cover,COVER_instance.weight_sep[i],COVER_instance.profit_sep,i);

				cover_val= cover_val + COVER_instance.magic_constant - COVER_instance.cover_size +1;

				for (int i = 0; i < COVER_instance.cover_size; i++)
				{
					COVER_instance.cover[i]=1-COVER_instance.cover[i];
				}


				if(cover_val > TOLL_SEP)
				{


					////////////////////////////////////////////////////////////////////////////////////////////////////////////
					double counter=0;
					for (int i = 0; i < MULTI_KP_instance.N_ITEMS; i++){if(COVER_instance.cover[i]>0.5){counter++;}}
					double viol=0;
					for (int i = 0; i < MULTI_KP_instance.N_ITEMS; i++){if(COVER_instance.cover[i]>0.5){viol+=(COVER_instance.point[i]-1);}}viol++;
					cout << "cover_val\t" << cover_val << "\tnumber_of_covers\t"<< COVER_instance.number_of_covers << "\t density \t" << counter/MULTI_KP_instance.N_ITEMS << "\t item \t" << counter << "\t viol \t" << viol << endl;
					viol_total+=viol;
					if(viol_max<viol){viol_max=viol;}
					if(viol_min>viol){viol_min=viol;}
					density_total+=counter/MULTI_KP_instance.N_ITEMS;
					if(density_max<counter/MULTI_KP_instance.N_ITEMS){density_max=counter/MULTI_KP_instance.N_ITEMS;}
					if(density_min>counter/MULTI_KP_instance.N_ITEMS){density_min=counter/MULTI_KP_instance.N_ITEMS;}
					item_total+=counter;
					if(item_max<counter){item_max=counter;}
					if(item_min>counter){item_min=counter;}
					////////////////////////////////////////////////////////////////////////////////////////////////////////////


					add_cover(&MULTI_KP_instance,COVER_instance.cover);

					COVER_instance.number_of_covers++;

					NOT_ADDED=false;
				}

			}

			if(NOT_ADDED)
			{
				break;
			}


		}
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////



		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		if(MULTI_KP_instance.algo==144)
		{

			bool NOT_ADDED=true;

			for(int i=0; i<MULTI_KP_instance.N_KNAPSACKS; i++)
			{

				cover_val=DP_kp01_advanced(COVER_instance.cover_size,COVER_instance.cap_sep[i], COVER_instance.profit_sep,COVER_instance.weight_sep[i],COVER_instance.cover,false,COVER_instance.order_DP[i]);

				make_maximal(&COVER_instance,COVER_instance.cover_size,COVER_instance.cap_sep[i],COVER_instance.cover,COVER_instance.weight_sep[i],COVER_instance.profit_sep,i);

				cover_val= cover_val + COVER_instance.magic_constant - COVER_instance.cover_size +1;

				for (int i = 0; i < COVER_instance.cover_size; i++)
				{
					COVER_instance.cover[i]=1-COVER_instance.cover[i];
				}


				if(cover_val > TOLL_SEP)
				{


					////////////////////////////////////////////////////////////////////////////////////////////////////////////
					double counter=0;
					for (int i = 0; i < MULTI_KP_instance.N_ITEMS; i++){if(COVER_instance.cover[i]>0.5){counter++;}}
					double viol=0;
					for (int i = 0; i < MULTI_KP_instance.N_ITEMS; i++){if(COVER_instance.cover[i]>0.5){viol+=(COVER_instance.point[i]-1);}}viol++;
					cout << "cover_val (EXTENDED)\t" << cover_val << "\tnumber_of_covers\t"<< COVER_instance.number_of_covers << "\t density \t" << counter/MULTI_KP_instance.N_ITEMS<< "\t item \t" << counter << "\tviol\t " << viol<< endl;
					viol_total+=viol;
					if(viol_max<viol){viol_max=viol;}
					if(viol_min>viol){viol_min=viol;}
					density_total+=counter/MULTI_KP_instance.N_ITEMS;
					if(density_max<counter/MULTI_KP_instance.N_ITEMS){density_max=counter/MULTI_KP_instance.N_ITEMS;}
					if(density_min>counter/MULTI_KP_instance.N_ITEMS){density_min=counter/MULTI_KP_instance.N_ITEMS;}
					item_total+=counter;
					if(item_max<counter){item_max=counter;}
					if(item_min>counter){item_min=counter;}
					////////////////////////////////////////////////////////////////////////////////////////////////////////////


					/////////////////////////////////////////////////////////////////////////////////////////////////////////
					int eureka=extend_cover(MULTI_KP_instance.N_ITEMS,MULTI_KP_instance.RHS_CAPACITIES[i],COVER_instance.cover,MULTI_KP_instance.MATRIX_WEIGHTS[i]);
					cout << "eureka\t" << eureka << endl;
					/////////////////////////////////////////////////////////////////////////////////////////////////////////

					add_cover_extended(&MULTI_KP_instance,COVER_instance.cover,counter);

					total_items_extended+=eureka;

					COVER_instance.number_of_covers++;
					NOT_ADDED=false;
				}

			}

			if(NOT_ADDED)
			{
				break;
			}


		}
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		if(MULTI_KP_instance.algo==7)
		{

			bool NOT_ADDED=true;

			for(int i=0; i<MULTI_KP_instance.N_KNAPSACKS; i++)
			{

				cover_val=solve_DP_lex(&COVER_instance,i);


				if(cover_val<1-TOLL_SEP)
				{


					////////////////////////////////////////////////////////////////////////////////////////////////////////////
					double counter=0;
					for (int i = 0; i < MULTI_KP_instance.N_ITEMS; i++){if(COVER_instance.cover[i]>0.5){counter++;}}
					double viol=0;
					for (int i = 0; i < MULTI_KP_instance.N_ITEMS; i++){if(COVER_instance.cover[i]>0.5){viol+=(COVER_instance.point[i]-1);}}viol++;
					cout << "cover_val\t" << cover_val << "\tnumber_of_covers\t"<< COVER_instance.number_of_covers << "\t density \t" << counter/MULTI_KP_instance.N_ITEMS << "\t item \t" << counter << "\t viol \t" << viol << endl;
					viol_total+=viol;
					if(viol_max<viol){viol_max=viol;}
					if(viol_min>viol){viol_min=viol;}
					density_total+=counter/MULTI_KP_instance.N_ITEMS;
					if(density_max<counter/MULTI_KP_instance.N_ITEMS){density_max=counter/MULTI_KP_instance.N_ITEMS;}
					if(density_min>counter/MULTI_KP_instance.N_ITEMS){density_min=counter/MULTI_KP_instance.N_ITEMS;}
					item_total+=counter;
					if(item_max<counter){item_max=counter;}
					if(item_min>counter){item_min=counter;}
					////////////////////////////////////////////////////////////////////////////////////////////////////////////


					add_cover(&MULTI_KP_instance,COVER_instance.cover);

					COVER_instance.number_of_covers++;

					NOT_ADDED=false;
				}

			}

			if(NOT_ADDED)
			{
				break;
			}


		}
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////



		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		if(MULTI_KP_instance.algo==77)
		{

			bool NOT_ADDED=true;

			for(int i=0; i<MULTI_KP_instance.N_KNAPSACKS; i++)
			{

				cover_val=solve_DP_lex(&COVER_instance,i);


				if(cover_val<1-TOLL_SEP)
				{


					////////////////////////////////////////////////////////////////////////////////////////////////////////////
					double counter=0;
					for (int i = 0; i < MULTI_KP_instance.N_ITEMS; i++){if(COVER_instance.cover[i]>0.5){counter++;}}
					double viol=0;
					for (int i = 0; i < MULTI_KP_instance.N_ITEMS; i++){if(COVER_instance.cover[i]>0.5){viol+=(COVER_instance.point[i]-1);}}viol++;
					cout << "cover_val (EXTENDED)\t" << cover_val << "\tnumber_of_covers\t"<< COVER_instance.number_of_covers << "\t density \t" << counter/MULTI_KP_instance.N_ITEMS<< "\t item \t" << counter << "\tviol\t " << viol<< endl;
					viol_total+=viol;
					if(viol_max<viol){viol_max=viol;}
					if(viol_min>viol){viol_min=viol;}
					density_total+=counter/MULTI_KP_instance.N_ITEMS;
					if(density_max<counter/MULTI_KP_instance.N_ITEMS){density_max=counter/MULTI_KP_instance.N_ITEMS;}
					if(density_min>counter/MULTI_KP_instance.N_ITEMS){density_min=counter/MULTI_KP_instance.N_ITEMS;}
					item_total+=counter;
					if(item_max<counter){item_max=counter;}
					if(item_min>counter){item_min=counter;}
					////////////////////////////////////////////////////////////////////////////////////////////////////////////


					/////////////////////////////////////////////////////////////////////////////////////////////////////////
					int eureka=extend_cover(MULTI_KP_instance.N_ITEMS,MULTI_KP_instance.RHS_CAPACITIES[i],COVER_instance.cover,MULTI_KP_instance.MATRIX_WEIGHTS[i]);
					cout << "eureka\t" << eureka << endl;
					/////////////////////////////////////////////////////////////////////////////////////////////////////////

					add_cover_extended(&MULTI_KP_instance,COVER_instance.cover,counter);

					total_items_extended+=eureka;

					COVER_instance.number_of_covers++;

					NOT_ADDED=false;
				}

			}

			if(NOT_ADDED)
			{
				break;
			}


		}
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////



	}



	//	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	// * writing the created ILP model on a file *
	//	MULTI_KP_instance.status=CPXwriteprob(MULTI_KP_instance.env_MULTI_KP,MULTI_KP_instance.lp_MULTI_KP,"FINAL.lp",NULL);
	//	if(MULTI_KP_instance.status!=0) {
	//		printf("error in CPXwriteprob\n");
	//		exit(-1);
	//	}
	//	/////////////////////////////////////////////////////////////////////////////////////////////////////////



	clock_t time_end=clock();
	double solution_time=(double)(time_end-time_start)/(double)CLOCKS_PER_SEC;


	multi_kp_free_cplex(&MULTI_KP_instance);

	for(int i=0; i<MULTI_KP_instance.N_KNAPSACKS; i++)
	{
		free_cplex_cover(&COVER_instance,i);
	}

	cout << "\n\nGENERATED\t" << COVER_instance.number_of_covers << "\tCOVERS" <<endl;
	cout << "TIME\t" << solution_time <<endl;
	cout << fixed << "LP_CURRENT\t" << LP_CURRENT << endl;

	cout << fixed << "Average density\t" << density_total/COVER_instance.number_of_covers << endl;
	cout << fixed << "Average max\t" << density_max<< endl;
	cout << fixed << "Average min\t" << density_min << endl;

	cout << fixed << "Average viol\t" << viol_total/COVER_instance.number_of_covers << endl;
	cout << fixed << "Average max\t" << viol_max<< endl;
	cout << fixed << "Average min\t" << viol_min << endl;


	cout << fixed << "Average items\t" << item_total/COVER_instance.number_of_covers << endl;
	cout << fixed << "Average items max\t" << item_max<< endl;
	cout << fixed << "Average items min\t" << item_min << endl;

	cout << fixed << "Average item_profit_frac\t" << total_item_profit_frac /total_iterations << endl;
	cout << fixed << "Average item_profit_zero\t" << total_item_profit_zero /total_iterations << endl;
	cout << fixed << "Average item_profit_uno\t" << total_item_profit_uno /total_iterations << endl;


	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	double max_i_avg=0;
	double cap_sep_avg=0;

	for(int i=0; i<MULTI_KP_instance.N_KNAPSACKS; i++)
	{
		max_i_avg+=COVER_instance.max_i[i];
		cap_sep_avg+=COVER_instance.cap_sep[i];
	}

	max_i_avg=max_i_avg/MULTI_KP_instance.N_KNAPSACKS;
	cap_sep_avg=cap_sep_avg/MULTI_KP_instance.N_KNAPSACKS;


	//	/////////////////////////////////////////////////////////////////////////////
	ofstream compact_file;
	compact_file.open("info_Exensive.txt", ios::app);
	compact_file << fixed



			<< MULTI_KP_instance.objval << "\t"
			<< MULTI_KP_instance.bestobjval << "\t"
			<< MULTI_KP_instance.solution_time << "\t"

			<< MULTI_KP_instance.stat << "\t"
			<< MULTI_KP_instance.nodecount << "\t"

			<<  MULTI_KP_instance.cplex_lp << "\t"
			<<  MULTI_KP_instance.solution_time_lp << "\t"



			<< COVER_instance.number_of_covers << "\t"
			<< solution_time << "\t"

			<<   LP_INIT << "\t"

			<<   LP_CURRENT << "\t"

			<<    max_i_avg << "\t"
			<<    cap_sep_avg << "\t"

			<< density_total/COVER_instance.number_of_covers << "\t"
			<< density_max<< "\t"
			<< density_min<< "\t"

			<< item_total/COVER_instance.number_of_covers << "\t"
			<< item_max<< "\t"
			<< item_min<< "\t"

			<< total_item_profit_frac /total_iterations << "\t"
			<< total_item_profit_zero /total_iterations << "\t"
			<< total_item_profit_uno /total_iterations << "\t"

			<< total_items_extended / COVER_instance.number_of_covers << "\t"

			<< viol_total/COVER_instance.number_of_covers << "\t"
			<< viol_max<< "\t"
			<< viol_min<< "\t"

			<< 			MULTI_KP_instance.n_conflicts << "\t"

			<< endl;
	compact_file.close();
	//	/////////////////////////////////////////////////////////////////////////////




	free(COVER_instance.env_sep);
	free(COVER_instance.lp_sep);


	COVER_DATA_free(&COVER_instance,MULTI_KP_instance.N_KNAPSACKS );

	printf("\nDONE!");

	free_data_prob_multi_kp(&MULTI_KP_instance);


	free(MULTI_KP_instance.instname);

	return 1;
}








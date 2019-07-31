#ifndef TEST_H
#define TEST_H
#include <omp.h>
#include "euler_field.h"
#include "flow_past_a_circular_cylinder.h"
void test04() {

	int Grid_X;
	int Grid_Y;
	double center_x;
	double center_y;
	double r;
	int number_of_node;

	double Re;
	double L;
	double u;
	int number_of_threads;
	int max_step;

	Grid_X = 800;
	Grid_Y = 400;

	center_x = double(Grid_X/4);
	center_y = double(Grid_Y/2);
	r = 50;
	number_of_node = 240;

	Re = 20;
	L = 2 * r;
	u = 0.1;
	number_of_threads = 32;
	max_step = 500001;


		flow_past_a_circular_cylinder flow1 = flow_past_a_circular_cylinder(Grid_X, Grid_Y,
			center_x, center_y, r, number_of_node,
			Re, L, u, number_of_threads, max_step);
		
		
		
		

		//flow1.multi_evolution();
		flow1.evolution();

		flow1.multi_evolution();
		

		//flow1.multi_direct_forcing_step();



		//flow1.Set_up_interreaction_data();
		
		
		cout << "end" << endl;





}


//Ô²ÖùÈÆÁ÷
//Re100-200

void test03() {
	int Grid_X;
	int Grid_Y;
	double center_x;
	double center_y;
	double r;
	int number_of_node;

	//double Re;
	double L;
	double u;
	int number_of_threads;
	int max_step;

	Grid_X = 500;
	Grid_Y = 80;

	center_x = 40;
	center_y = 125;
	r = 10;
	number_of_node = 500;

	//Re = 20;
	L = 2 * r;
	u = 0.1;
	number_of_threads = 16;
	max_step = 150001;


	for (int Re = 50; Re < 201; Re = Re + 50)
	{
		flow_past_a_circular_cylinder flow1 = flow_past_a_circular_cylinder(Grid_X, Grid_Y,
			center_x, center_y, r, number_of_node,
			Re, L, u, number_of_threads, max_step);
		flow1.evolution();
		cout << "end" << endl;
	}



}

//Ô²ÖùÈÆÁ÷
//Re1-70
void test02() {

	int Grid_X = 20;
	int Grid_Y = 20;
	double center_x=double(Grid_X/2);
	double center_y=double(Grid_Y/2);
	double r=5;
	int number_of_node = 200;

	double Re=0;
	double L= 2*r;
	double u=0.1;
	int number_of_threads=4;
	int max_step=101;




	for (int i = 5; i < 70; i=i+5)
	{
		Re = i;
		flow_past_a_circular_cylinder flow1 = flow_past_a_circular_cylinder(Grid_X, Grid_Y,
			center_x, center_y, r, number_of_node,
			Re, L, u, number_of_threads, max_step);
		flow1.evolution();
		cout << "end" << endl;
	}





}


void test01() {
	omp_set_num_threads(16);
	std::cout << "threads =" << omp_get_num_threads() << endl;
	int grid_x;
	int grid_y;
	double Re;
	double u;
	double L;
	int max_step;

	Re = 50;
	grid_x = 10;
	grid_y = 10;
	L = double(2 * grid_y);

	u = 0.1;
	max_step = 10000;
	euler_field Poiseuille = euler_field(grid_x, grid_y, Re, L, u, max_step);
	Poiseuille.evolution();



}






#endif // TEST_H

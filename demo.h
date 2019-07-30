#include "flow_past_a_circular_cylinder.h"

/*lpl*/
void demo1() {

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

	int NF = 20;
	Grid_X =  100;
	Grid_Y =  100;

	center_x = double(Grid_X / 4);
	center_y = double(Grid_Y / 2);
	r = 20;

	number_of_node = 30;

	Re = 20;
	L = 2 * r;
	u = 0.1;
	number_of_threads = 24;
	max_step = 3001;


	flow_past_a_circular_cylinder flow1 = flow_past_a_circular_cylinder(Grid_X, Grid_Y,
		center_x, center_y, r, number_of_node,
		Re, L, u, number_of_threads, max_step ,NF);
	
	//flow1.mluti_method_evolution();


	flow1.evolution();



	//flow1.clear_NF_Parameter();


	//flow1.mluti_test_step_a();
	//flow1.mluti_test_step_b();

	//


	//flow1.mluti_test_step_c(1);
	//flow1.mluti_test_step_d(1);
	//flow1.mluti_test_step_e(1);
	//flow1.mluti_test_step_f(1);



	//flow1.mluti_test_step_c(2);
	//flow1.mluti_test_step_d(2);
	//flow1.mluti_test_step_e(2);
	//flow1.mluti_test_step_f(2);




	//flow1.mluti_test_step_c(3);
	//flow1.mluti_test_step_d(3);
	//flow1.mluti_test_step_e(3);
	//flow1.mluti_test_step_f(3);


	//flow1.mluti_test_step_c(4);
	//flow1.mluti_test_step_d(4);
	//flow1.mluti_test_step_e(4);
	//flow1.mluti_test_step_f(4);

	//




	//flow1.mluti_method_step_a();
	//flow1.mluti_method_step_b();

	




	//for (int m = 1; m < 5; m++)
	//{
	//	flow1.mluti_test_step_c(m);

	//	flow1.mluti_test_step_d(m);

	//	flow1.mluti_method_step_e(m);

	//	flow1.mluti_test_step_f(m);
	//

	//}


	//for (int i = 0; i < 500; i++)
	//{
	//	flow1.mluti_method_first_forcing_step();
	//	//flow1.mluti_test_step_g();


	//	flow1.Collision_step();
	//	flow1.Secord_forcing_step();
	//	flow1.Stream_step();

	//	flow1.Macroscopic();

	//}



	//Êä³ö
	//flow1.Output_mluti_test();


	//system("pause");
	cout << "end" << endl;







}






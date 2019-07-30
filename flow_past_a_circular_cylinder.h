#ifndef FLOW_PAST_A_CIRCULAR_CYLINDER_H
#define FLOW_PAST_A_CIRCULAR_CYLINDER_H
#include "Lagrange_field.h"
#include "euler_field.h"
#include "Interreaction.h"
#include "Function.h"


class flow_past_a_circular_cylinder
{
private:
	//euler node
	int X_Max;
	int Y_Max;

	//lagrange node
	double _r;
	double _center_x;
	double _center_y;
	int _Number_of_node;





	//Parameter
	double _u;
	double _Re;
	double _L;
	//double _tau;
	int _number_of_threads;
	int _max_step;

	//lagrange_area

	int _lagrange_area_x1;
	int _lagrange_area_x2;

	int _lagrange_area_y1;
	int _lagrange_area_y2;



	Lagrange_field* Solid;
	euler_field* Fluid;

	Interreaction* _interreaction;


public:
	flow_past_a_circular_cylinder();

	flow_past_a_circular_cylinder(int Grid_X, int Grid_Y,
		double center_x, double center_y, double r, int number_of_node, 
		double Re, double L, double u, int number_of_threads,int max_step);

	//euler node
	void Set_X_Max(int x);
	void Set_Y_Max(int y);

	//lagrange node
	void Set_number_of_node(int N);
	void Set_center_x(double x);
	void Set_center_y(double y);
	void Set_r(double r);

	//Parameter
	void Set_Re(double Re);
	void Set_L(double L);
	void Set_u(double u);
	void Set_number_of_threads(int N);
	//void Set_tau();
	void Set_max_step(int Step);

	//lagrange_area
	void Set_lagrange_area();


	void Set_Solid();
	void Set_Fuild();

	void Set_Interreaction();

	double Get_X_Max();
	double Get_Y_Max();

	double Get_center_x();
	double Get_center_y();
	double Get_r();

	double Get_Re();
	double Get_L();
	double Get_number_of_threads();
	double Get_u();

	Lagrange_field* Get_Solid();
	euler_field* Get_Fuild();




	//直接力法步骤

	//step a
	void u_ij_step();
	//step b
	void Unforce_velocity_interpolation();
	//step c
	void Boundary_force_evaluation_on_Xb();
	//step d
	void Force_distribution_on_Xij();
	//step e
	void Update_of_velocity_on_xij();



	//力项lbm步骤

	void First_forcing_step();

	void Collision_step();

	void Secord_forcing_step();

	void Stream_step();

	void Macroscopic();

	void evolution();

	//set_interreaction_data

	void Set_up_interreaction_data();

	//multi-direct-forcing

	
	void multi_evolution();

	void multi_direct_forcing_step();

	void euler_node_to_interreaction();

	void lagrange_node_to_interreation();

	void interreaction_to_lagrange_node();

	void interreaction_to_euler_node();


	void Output_Parameter();
	void Output_Data();
	void Output_Solid(double value,int step);
	void Output_Fluid(double value, int step);





	~flow_past_a_circular_cylinder();
};

#endif // _DEBUG
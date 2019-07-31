#ifndef INTERREACTION
#define INTERREACTION
#include "Temporary_Euler_node.h"
#include "Temporary_Lagrange_Node.h"

class Interreaction
{
public:
	int _nf;

	int _euler_node_x_max;
	int _euler_node_y_max;

	int _lagrange_node_number;

	int _euler_node_number;



	vector<vector<bool> > _Is_Compute;

	vector<Temporary_Euler_node*>  _Euler_data;
	
	vector<Temporary_Lagrange_Node*> _lagrange_data;


public:
	Interreaction(int x_max ,int y_max,int lagrange_node_number);

	Interreaction();
	
	void Set_euler_node_number();
	void Set_lagrange_node_number();

	void init_euler_data();
	void init_lagrange_data();
	
	void Set_valid_cal();



	void Output_euler_data();
	void Output_laghrange_data();


	void Evolution_Step();

	void Multi_step_a();
	
	void Multi_step_b();

	void Multi_step_c(int m);

	void Multi_step_d(int m);

	void Multi_step_e(int m);

	void Multi_step_f(int m);
	
	void Calculate_result();

	void output_euler_node(double value,int step);


	~Interreaction();



	//void Interreaction_of_Solid_and_Fluid();




};

#endif // !INTERREACTION

#ifndef TEMPORARY_LAGRANGE_NODE
#define TEMPORARY_LAGRANGE_NODE
#include "euler_field.h"
#include "Lagrange_field.h"
#include "Temporary_Euler_node.h"

class Temporary_Lagrange_Node
{
public:

	int _nf;

	int _node_number;

	double _delta_sb;

	double _grid_x;
	double _grid_y;
	//Lagrange_Point* _temporary_lagrange_node;

	double _U_b_x;
	double _U_b_y;

	vector<double> _u_b_x;
	vector<double> _u_b_y;



	vector<double> _F_b_x;
	vector<double> _F_b_y;

	double _result_F_b_x;
	double _result_F_b_y;

	//double _result_ub_x;
	//double _result_ub_y;


	vector<int> _valid_euler_node;
	//vector<Temporary_Euler_node*> _valid_euler_node;
	//vector<int> valid_y;


public:
	Temporary_Lagrange_Node();

	//Temporary_Lagrange_Node(double grid_x,double grid_y,
	//	double Ub_x,double Ub_y);
	
	//Temporary_Lagrange_Node(double grid_x, double grid_y);
	
	Temporary_Lagrange_Node(int node_number,int nf);
	
	
	//void Set_Ub(double u,double v);

	void Set_delta_sb(double value);

	void Set_nf(int nf);

	void Set_node_number(int a);

	void Set_grid_xy(double x,double y);

	void Set_U_b(double x, double y);

	void Set_u_b(double u, double v ,int m);

	void Set_F_b(double fx, double fy, int m);

	void Add_valid_node(int x);

	void Output_infomation_in_screen();

	double Get_position_x();
	double Get_position_y();




	int Get_valid_euler_node_size();
	int Get_valid_euler_node(int node);

	double Get_U_b_x();
	double Get_U_b_y();

	double Get_F_b_x(int a);
	double Get_F_b_y(int a);

	double Get_u_b_x(int a);
	double Get_u_b_y(int a);

	double Get_delta_sb();

	double Get_result_F_b_x();
	double Get_result_F_b_y();

	int Get_node_number();

	~Temporary_Lagrange_Node();
};

#endif // !Temporary_Lagrange_Node

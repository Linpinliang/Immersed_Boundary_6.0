#ifndef TEMPORARY_EULER_NODE
#define TEMPORARY_EULER_NODE
#include "euler_field.h"
#include "Lagrange_field.h"
#include "Temporary_Lagrange_Node.h"

class Temporary_Euler_node
{
public:
	bool _calculate;
	int _nf = 20;

	int _number_x;
	int _number_y;

	double _grid_x;
	double _grid_y;

	

	//Euler_Point* _temporary_euler_node;

	vector<double> _u_ij_x;
	vector<double> _u_ij_y;

	vector<double> _F_ij_x;
	vector<double> _F_ij_y;

	double _result_F_ij_x;
	double _result_F_ij_y;

	double _result_u_ij_x;
	double _result_u_ij_y;

	vector<int> _valid_lagrange_node;

public:

	Temporary_Euler_node();

	Temporary_Euler_node(int number_x,int number_y,int nf);

	void Set_nf(int nf);

	
	void Set_u_ij_m(double u, double v, int m);

	void Set_F_ij_M(double fx, double fy, int m);

	void Set_Position(double x ,double y);

	void add_valid_node(int node);


	double Get_Position_x();
	double Get_Position_y();

	double Get_number_x();
	double Get_number_y();

	double Get_u_ij_x(int a);
	double Get_u_ij_y(int a);

	double Get_F_ij_x(int a);
	double Get_F_ij_y(int a);



	double Get_result_F_ij_x();
	double Get_result_F_ij_y();

	double Get_result_u_ij_x();
	double Get_result_u_ij_y();

	int Get_valid_node_size();
	int Get_valid_node_number(int node);


	void output_data();



	~Temporary_Euler_node();
};


#endif // !TEMPORARY_EULER_NODE

#include "Temporary_Euler_node.h"



Temporary_Euler_node::Temporary_Euler_node()
{
}

Temporary_Euler_node::Temporary_Euler_node(int number_x, int number_y, int nf)
{
	Set_nf(nf);
	_number_x = number_x;
	_number_y = number_y;


}




void Temporary_Euler_node::Set_nf(int nf)
{
	_nf = nf;

	//_u_ij_x.reserve(nf);
	//_u_ij_y.reserve(nf);

	_u_ij_x = vector<double>(nf);
	_u_ij_y = vector<double>(nf);

	_F_ij_x = vector<double>(nf);
	_F_ij_y = vector<double>(nf);


}

void Temporary_Euler_node::Set_u_ij_m(double u, double v, int m)
{


	_u_ij_x[m] = u;
	_u_ij_y[m] = v;


}

void Temporary_Euler_node::Set_F_ij_M(double fx, double fy, int m)
{
	_F_ij_x[m] = fx;
	_F_ij_y[m] = fy;


}

void Temporary_Euler_node::Set_Position(double x, double y)
{
	_grid_x = x;
	_grid_y = y;

}

void Temporary_Euler_node::add_valid_node(int node)
{
	_valid_lagrange_node.push_back(node);
}



double Temporary_Euler_node::Get_Position_x()
{
	return _grid_x;
}

double Temporary_Euler_node::Get_Position_y()
{
	return _grid_y;
}

double Temporary_Euler_node::Get_number_x()
{
	return _number_x;
}

double Temporary_Euler_node::Get_number_y()
{
	return _number_y;
}

double Temporary_Euler_node::Get_u_ij_x(int a)
{
	return _u_ij_x[a];
}

double Temporary_Euler_node::Get_u_ij_y(int a)
{
	return _u_ij_y[a];
}

double Temporary_Euler_node::Get_F_ij_x(int a)
{
	return _F_ij_x[a];
}

double Temporary_Euler_node::Get_F_ij_y(int a)
{
	return _F_ij_y[a];
}



double Temporary_Euler_node::Get_result_u_ij_x()
{

	_result_u_ij_x = _u_ij_x[_nf-1];


	return _result_u_ij_x;
}

double Temporary_Euler_node::Get_result_u_ij_y()
{

	_result_u_ij_y = _u_ij_y[_nf - 1];

	return _result_u_ij_y;
}

int Temporary_Euler_node::Get_valid_node_size()
{
	return _valid_lagrange_node.size();
}

int Temporary_Euler_node::Get_valid_node_number(int node)
{
	return _valid_lagrange_node[node];
}

void Temporary_Euler_node::output_data()
{
	cout << _number_x<<"\t"<<_number_y << endl;
	cout << _grid_x << "\t" << _grid_y << endl;


	for (int i = 0; i < _u_ij_x.size(); i++)
	{
		cout << _u_ij_x[i] << "\t";
	}

	cout << endl;
	//cout << _u_ij_x[0];

	for (int i = 0; i < _u_ij_x.size(); i++)
	{
		cout << _u_ij_y[i] << "\t";
	}

	cout << endl;

	for (int i = 0; i < _u_ij_x.size(); i++)
	{
		cout << _F_ij_x[i] << "\t";
	}


	cout << endl;
	for (int i = 0; i < _u_ij_x.size(); i++)
	{
		cout << _F_ij_y[i] << "\t";
	}


	cout << endl;

	//cout << "valid lagrange node:" << _valid_lagrange_node.size();

	//cout << endl;

	cout << _result_u_ij_x << "\t" << _result_u_ij_y;

	cout << endl; 
	
	cout << _result_F_ij_x << "\t" << _result_F_ij_y;
	
	
	cout << endl;
	


}



double Temporary_Euler_node::Get_result_F_ij_x()
{
	double a = 0;

	for (int i = 0; i < _F_ij_x.size(); i++)
	{
		a += _F_ij_x[i];

	}

	_result_F_ij_x = a;

	return a;
}

double Temporary_Euler_node::Get_result_F_ij_y()
{
	double a = 0;

	for (int i = 0; i < _F_ij_y.size(); i++)
	{
		a += _F_ij_y[i];

	}


	_result_F_ij_y = a;

	return a;
}

Temporary_Euler_node::~Temporary_Euler_node()
{
}

#include "Temporary_Lagrange_Node.h"



Temporary_Lagrange_Node::Temporary_Lagrange_Node()
{
}



Temporary_Lagrange_Node::Temporary_Lagrange_Node(int node_number, int nf)
{
	Set_nf(nf);
	_node_number = node_number;


}


void Temporary_Lagrange_Node::Set_delta_sb(double value)
{
	_delta_sb = value;
}

void Temporary_Lagrange_Node::Set_nf(int nf)
{
	_nf = nf;

	//vector<double> _u_b_x;
	//vector<double> _u_b_y;

	//vector<double> _F_b_x;
	//vector<double> _F_b_y;

	_u_b_x = vector<double>(nf);
	_u_b_y = vector<double>(nf);
	_F_b_x = vector<double>(nf);
	_F_b_y = vector<double>(nf);

}

void Temporary_Lagrange_Node::Set_node_number(int a)
{
	_node_number = a;
}

void Temporary_Lagrange_Node::Set_grid_xy(double x, double y)
{
	_grid_x = x;
	_grid_y = y;

}

void Temporary_Lagrange_Node::Set_U_b(double x, double y)
{
	_U_b_x = x;
	_U_b_y = y;

}

void Temporary_Lagrange_Node::Set_u_b(double u, double v, int m)
{
	_u_b_x[m] = u;
	_u_b_y[m] = v;


}

void Temporary_Lagrange_Node::Set_F_b(double fx, double fy, int m)
{
	_F_b_x[m] = fx;
	_F_b_y[m] = fy;

}

void Temporary_Lagrange_Node::Add_valid_node(int x)
{
	_valid_euler_node.push_back(x);
}

void Temporary_Lagrange_Node::Output_infomation_in_screen()
{
	cout << _node_number<<endl;

	cout << _delta_sb << endl;

	cout << _grid_x << "\t" << _grid_y << endl;
	cout << _U_b_x << "\t" << _U_b_y << endl;

	for (int i = 0; i < _u_b_x.size(); i++)
	{
		cout << _u_b_x[i] << "\t";
	}

	cout << endl;

	for (int i = 0; i < _u_b_x.size(); i++)
	{
		cout << _u_b_y[i] << "\t";
	}


	cout << endl;


	for (int i = 0; i < _u_b_x.size(); i++)
	{
		cout << _F_b_x[i] << "\t";
	}


	cout << endl;


	for (int i = 0; i < _u_b_x.size(); i++)
	{
		cout << _F_b_y[i] << "\t";
	}

	/*
	cout << endl;

	cout << "nodesize():" << _valid_euler_node.size() << endl;

	for (int i = 0; i < _valid_euler_node.size(); i++)
	{
		cout << _valid_euler_node[i]<<"\t";
	}

	*/



	cout << endl;


	cout << _result_F_b_x << "\t" << _result_F_b_y;

	cout << endl;



}

double Temporary_Lagrange_Node::Get_position_x()
{
	return _grid_x;
}

double Temporary_Lagrange_Node::Get_position_y()
{
	return _grid_y;
}

int Temporary_Lagrange_Node::Get_valid_euler_node_size()
{
	return _valid_euler_node.size();
}

int Temporary_Lagrange_Node::Get_valid_euler_node(int node)
{
	return _valid_euler_node[node];
}

double Temporary_Lagrange_Node::Get_U_b_x()
{
	return _U_b_x;
}

double Temporary_Lagrange_Node::Get_U_b_y()
{
	return _U_b_y;
}

double Temporary_Lagrange_Node::Get_F_b_x(int a)
{
	return _F_b_x[a];
}

double Temporary_Lagrange_Node::Get_F_b_y(int a)
{
	return _F_b_y[a];
}

double Temporary_Lagrange_Node::Get_u_b_x(int a)
{
	return _u_b_x[a];
}

double Temporary_Lagrange_Node::Get_u_b_y(int a)
{
	return _u_b_y[a];
}

double Temporary_Lagrange_Node::Get_delta_sb()
{
	return _delta_sb;
}

double Temporary_Lagrange_Node::Get_result_F_b_x()
{
	_result_F_b_x = 0;
	
	
	for (int i = 0; i < _F_b_x.size(); i++)
	{
		_result_F_b_x += _F_b_x[i];
	}



	return _result_F_b_x;
}

double Temporary_Lagrange_Node::Get_result_F_b_y()
{
	_result_F_b_y = 0;

	for (int i = 0; i < _F_b_y.size(); i++)
	{
		_result_F_b_y += _F_b_y[i];


	}



	return  _result_F_b_y;
}

int Temporary_Lagrange_Node::Get_node_number()
{
	return _node_number;
}

Temporary_Lagrange_Node::~Temporary_Lagrange_Node()
{
}

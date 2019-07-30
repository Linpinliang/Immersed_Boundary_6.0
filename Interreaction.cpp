#include "Interreaction.h"



Interreaction::Interreaction(int x_max, int y_max, int lagrange_node_number)
{
	_euler_node_x_max = x_max;
	_euler_node_y_max = y_max;

	_lagrange_node_number = lagrange_node_number;

	//int _nf = 5;


	//_Euler_data = vector<vector<Temporary_Euler_node*> >(_euler_node_x_max + 1, vector<Temporary_Euler_node*>( _euler_node_y_max + 1));

	//_lagrange_data = vector<Temporary_Lagrange_Node*>(lagrange_node_number);
	
	_Is_Compute = vector<vector<bool> >(_euler_node_x_max + 1, vector<bool>(_euler_node_y_max + 1));



}

Interreaction::Interreaction()
{
}

void Interreaction::Set_euler_node_number()
{
	_euler_node_number = 0;


	for (int i = 0; i < _Is_Compute.size(); i++)
	{
		for (int j = 0; j < _Is_Compute[0].size(); j++)
		{
			if (_Is_Compute[i][j] == true) {
				_euler_node_number += 1;			

				//cout << i <<","<< j << endl;
				//cout << _euler_node_number << endl;
			}
		}
	}

	//cout<<"1euler_node_number:"<<_euler_node_number<<endl;

	_Euler_data = vector<Temporary_Euler_node*>(_euler_node_number);

	//cout << "2euler_node_number:" << _euler_node_number << endl;




}

void Interreaction::Set_lagrange_node_number()
{
	_lagrange_data = vector<Temporary_Lagrange_Node*>(_lagrange_node_number);
}

void Interreaction::init_euler_data()
{
	int counter = 0;

	for (int i = 0; i < _Is_Compute.size(); i++)
	{
		for (int j = 0; j < _Is_Compute[0].size(); j++)
		{
			if (_Is_Compute[i][j] == true) {

				Temporary_Euler_node* p1 = new Temporary_Euler_node(i, j, _nf);

				_Euler_data[counter] = p1;

				counter += 1;

			}
		}
	}



	//cout <<"euler_node_size()" <<_Euler_data.size()<<endl;

}

void Interreaction::init_lagrange_data()
{
	for (int i = 0; i < _lagrange_node_number; i++)
	{
		Temporary_Lagrange_Node* p1 = new Temporary_Lagrange_Node(i,_nf);

		_lagrange_data[i] = p1;
	}

	//cout << "lagrange_init()" << endl;
	//cout << "_lagrange_data.size():"<<_lagrange_data.size() << endl;

}

void Interreaction::Set_valid_cal()
{
	for (int i = 0; i < _Euler_data.size(); i++)
	{
		double ex = _Euler_data[i]->Get_Position_x();
		double ey = _Euler_data[i]->Get_Position_y();

		for (int j = 0; j < _lagrange_data.size(); j++)
		{
			double lx = _lagrange_data[j]->Get_position_x();
			double ly = _lagrange_data[j]->Get_position_y();

			if (D_function(ex, ey, lx, ly) > 0)
			{
				//cout << ex << "\t" << ey << "\t" << lx << "\t" << ly << endl;
				_Euler_data[i]->add_valid_node(j);
				_lagrange_data[j]->Add_valid_node(i);


			}



		}



	}

}

void Interreaction::Output_euler_data()
{	
	for (int i = 0; i < _Euler_data.size(); i++)
	{
		_Euler_data[i]->output_data();

	}



}

void Interreaction::Output_laghrange_data()
{
	for (int i = 0; i < _lagrange_data.size(); i++)
	{
		_lagrange_data[i]->Output_infomation_in_screen();
	}

}

void Interreaction::Evolution_Step()
{
	Multi_step_a();

	Multi_step_b();

	for (int m = 1; m < _nf; m++)
	{
		Multi_step_c(m);

		Multi_step_d(m);

		Multi_step_e(m);

		Multi_step_f(m);
	}


	Calculate_result();


	/*

	Multi_step_c(1);

	Multi_step_d(1);

	Multi_step_e(1);

	Multi_step_f(1);

	Multi_step_c(2);

	Multi_step_d(2);

	Multi_step_e(2);

	Multi_step_f(2);
	*/



}

void Interreaction::Multi_step_a()
{






}

void Interreaction::Multi_step_b()
{
	/*
	for (int node = 0; node < _lagrange_data.size(); node++)
	{

		double ub = 0;
		double vb = 0;

		for (int e_node_number = 0; e_node_number < _Euler_data.size(); e_node_number++)
		{
			double ux = _Euler_data[e_node_number]->_u_ij_x[0];
			double uy = _Euler_data[e_node_number]->_u_ij_y[0];

			double ex = _Euler_data[e_node_number]->Get_Position_x();
			double ey = _Euler_data[e_node_number]->Get_Position_y();

			double lx = _lagrange_data[node]->Get_position_x();
			double ly = _lagrange_data[node]->Get_position_y();


			//double D = Distance(ex, ey, lx, ly);
			double D = D_function(ex, ey, lx, ly);
			
			double h = 1;

			//
			
			if (D <= 2) {
				cout << "ex:" << ex << endl;
				cout << "ey:" << ey << endl;

				cout << "lx:" << lx << endl;
				cout << "ly:" << ly << endl;
			
				cout << "-----------------" <<endl;
			}

			
			//


			ub += ux * D * h * h;
			vb += uy * D * h * h;

		}
		cout << "interreaction_ub_vb:"<< ub << vb <<endl;


		_lagrange_data[node]->Set_u_b(ub, vb, 0);

	}

	*/

	#pragma omp parallel for
	for (int node = 0; node < _lagrange_data.size(); node++)
	{

		double ub = 0;
		double vb = 0;

		for (int i = 0; i < _lagrange_data[node]->Get_valid_euler_node_size(); i++)
		{
			int e_node_number = _lagrange_data[node]->Get_valid_euler_node(i);

			double ux = _Euler_data[e_node_number]->_u_ij_x[0];
			double uy = _Euler_data[e_node_number]->_u_ij_y[0];

			double ex = _Euler_data[e_node_number]->Get_Position_x();
			double ey = _Euler_data[e_node_number]->Get_Position_y();

			double lx = _lagrange_data[node]->Get_position_x();
			double ly = _lagrange_data[node]->Get_position_y();


			//double D = Distance(ex, ey, lx, ly);
			double D = D_function(ex, ey, lx, ly);

			double h = 1;



			ub += ux * D * h * h;
			vb += uy * D * h * h;

		}
		//cout << "interreaction_ub_vb:" << ub << vb << endl;


		_lagrange_data[node]->Set_u_b(ub, vb, 0);

	}


}

void Interreaction::Multi_step_c(int m)
{
	#pragma omp parallel for
	for (int i = 0; i < _lagrange_data.size(); i++)
	{
		double Ubx = _lagrange_data[i]->Get_U_b_x();
		double Uby = _lagrange_data[i]->Get_U_b_y();

		double ubx = _lagrange_data[i]->Get_u_b_x(m-1);
		double uby = _lagrange_data[i]->Get_u_b_y(m-1);

		double delta_t = 1;

		double rho = 1;

		double Fbx_m = 2 * rho / delta_t * (Ubx - ubx);
		double Fby_m = 2 * rho / delta_t * (Uby - uby);

		_lagrange_data[i]->Set_F_b(Fbx_m, Fby_m, m);
		//_lagrange_data[i]->Set_F_b(1, 2, m);

	}





}

void Interreaction::Multi_step_d(int m)
{
#pragma omp parallel for
	for (int i = 0; i < _Euler_data.size(); i++)
	{
		double ex = _Euler_data[i]->Get_Position_x();
		double ey = _Euler_data[i]->Get_Position_y();

		double F_ij_x = 0;
		double F_ij_y = 0;


		for (int j = 0; j < _Euler_data[i]->Get_valid_node_size(); j++)
		{
			int l_node = _Euler_data[i]->Get_valid_node_number(j);


			double lx = _lagrange_data[l_node]->Get_position_x();
			double ly = _lagrange_data[l_node]->Get_position_y();


			double D = D_function(ex, ey, lx, ly);
			double delta_sb = _lagrange_data[l_node]->Get_delta_sb();

			double Fbx = _lagrange_data[l_node]->Get_F_b_x(m);
			double Fby = _lagrange_data[l_node]->Get_F_b_y(m);

			F_ij_x += Fbx * D * delta_sb;
			F_ij_y += Fby * D * delta_sb;

		}


		_Euler_data[i]->Set_F_ij_M(F_ij_x, F_ij_y, m);



	}





}

void Interreaction::Multi_step_e(int m)
{
#pragma omp parallel for
	for (int i = 0; i < _Euler_data.size(); i++)
	{
		double delta_t = 1;
		double rho = 1;


		double u_before_x = _Euler_data[i]->Get_u_ij_x(m - 1);
		double u_before_y = _Euler_data[i]->Get_u_ij_y(m - 1);

		double F_ij_x = _Euler_data[i]->Get_F_ij_x(m);
		double F_ij_y = _Euler_data[i]->Get_F_ij_y(m);

		double ux = u_before_x + delta_t * F_ij_x / (2 * rho);
		double uy = u_before_y + delta_t * F_ij_y / (2 * rho);

		_Euler_data[i]->Set_u_ij_m(ux, uy, m);
	}




}





void Interreaction::Multi_step_f(int m)
{
	#pragma omp parallel for
	for (int i = 0; i < _lagrange_data.size(); i++)
	{
		double ubx = 0;
		double uby = 0;

		double lx = _lagrange_data[i]->Get_position_x();
		double ly = _lagrange_data[i]->Get_position_y();


		for (int j = 0; j < _lagrange_data[i]->Get_valid_euler_node_size(); j++)
		{
			int e_valid = _lagrange_data[i]->Get_valid_euler_node(j);

			double ex = _Euler_data[e_valid]->Get_Position_x();
			double ey = _Euler_data[e_valid]->Get_Position_y();

			double D = D_function(lx, ly, ex, ey);

			double h = 1;

			double u_ij_x_m = _Euler_data[e_valid]->Get_u_ij_x(m);
			double u_ij_y_m = _Euler_data[e_valid]->Get_u_ij_y(m);

			ubx += u_ij_x_m * D * h * h;
			uby += u_ij_y_m * D * h * h;



		}

		_lagrange_data[i]->Set_u_b(ubx, uby, m);



	}




}

void Interreaction::Calculate_result()
{

	
#pragma omp parallel for
	for (int i = 0; i < _Euler_data.size(); i++)
	{
		_Euler_data[i]->Get_result_F_ij_x();
		_Euler_data[i]->Get_result_F_ij_y();
		_Euler_data[i]->Get_result_u_ij_x();
		_Euler_data[i]->Get_result_u_ij_y();

	}



#pragma omp parallel for
	for (int i = 0; i <_lagrange_data.size(); i++)
	{
		_lagrange_data[i]->Get_result_F_b_x();
		_lagrange_data[i]->Get_result_F_b_y();
	}



	

}

void Interreaction::output_euler_node(double value, int step)
{
	ostringstream name2;
	name2 << "Re=" << value << "Interreaction_euler_node_" << step << ".plt";
	ofstream out(name2.str().c_str());

	void output_euler_node(double value, int step);
	
	for (int i = 0; i < _Euler_data.size(); i++)
	{
		out << _Euler_data[i]->Get_number_x() << "\t" << _Euler_data[i]->Get_number_y() << endl;
		out << _Euler_data[i]->_grid_x << "\t" << _Euler_data[i]->_grid_y << endl;


		for (int j = 0; j < _Euler_data[i]->_u_ij_x.size(); j++)
		{
			out << _Euler_data[i]->_u_ij_x[j] << "\t";
		}

		out << endl;
		//cout << _u_ij_x[0];

		for (int j = 0; j < _Euler_data[i]->_u_ij_x.size(); j++)
		{
			out << _Euler_data[i]->_u_ij_y[j] << "\t";
		}

		out << endl;

		for (int j = 0; j < _Euler_data[i]->_u_ij_x.size(); j++)
		{
			out << _Euler_data[i]->Get_F_ij_x(j) << "\t";
			//out << "aaa";
		}


		out << endl;
		for (int j = 0; j < _Euler_data[i]->_u_ij_x.size(); j++)
		{
			out << _Euler_data[i]->Get_F_ij_y(j) << "\t";
			//out << "aaa";
		}


		out << endl;

		//cout << "valid lagrange node:" << _valid_lagrange_node.size();

		//cout << endl;

		out << _Euler_data[i]->_result_u_ij_x << "\t" << _Euler_data[i]->_result_u_ij_y;

		out << endl;

		out << _Euler_data[i]->_result_F_ij_x << "\t" << _Euler_data[i]->_result_F_ij_y;


		out << endl;




	}



}





Interreaction::~Interreaction()
{
	//delete _Euler_data;
	//vector<Temporary_Euler_node*>  _Euler_data;

	//vector<Temporary_Lagrange_Node*> _lagrange_data;


	for (int i = 0; i < _Euler_data.size(); i++)
	{
		delete _Euler_data[i];

	}

	for (int i = 0; i < _lagrange_data.size(); i++)
	{
		delete _lagrange_data[i];
	}


}



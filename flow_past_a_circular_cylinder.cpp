#include "flow_past_a_circular_cylinder.h"




flow_past_a_circular_cylinder::flow_past_a_circular_cylinder()
{
}

flow_past_a_circular_cylinder::flow_past_a_circular_cylinder(int Grid_X, int Grid_Y, 
	double center_x, double center_y, double r, int number_of_node,
	double Re, double L, double u,int number_of_threads,int max_step)
{
	Set_X_Max(Grid_X);
	Set_Y_Max(Grid_Y);

	Set_number_of_node(number_of_node);
	Set_center_x(center_x);
	Set_center_y(center_y);
	Set_r(r);

	Set_Re(Re);
	Set_u(u);
	Set_L(L);
	Set_number_of_threads(number_of_threads);
	Set_max_step(max_step);


	Set_lagrange_area();


	Set_Solid();
	Set_Fuild();

}

void flow_past_a_circular_cylinder::Set_X_Max(int x)
{
	X_Max = x;
}

void flow_past_a_circular_cylinder::Set_Y_Max(int y)
{
	Y_Max = y;
}

void flow_past_a_circular_cylinder::Set_number_of_node(int N)
{
	_Number_of_node = N;
}

void flow_past_a_circular_cylinder::Set_center_x(double x)
{
	_center_x = x;
}

void flow_past_a_circular_cylinder::Set_center_y(double y)
{
	_center_y = y;
}

void flow_past_a_circular_cylinder::Set_r(double r)
{
	_r = r;
}

void flow_past_a_circular_cylinder::Set_Re(double Re)
{
	_Re = Re;
}

void flow_past_a_circular_cylinder::Set_L(double L)
{
	_L = L;

}

void flow_past_a_circular_cylinder::Set_u(double u)
{
	_u = u;
}

void flow_past_a_circular_cylinder::Set_number_of_threads(int N)
{
	_number_of_threads = N;

	omp_set_num_threads(_number_of_threads);

}

void flow_past_a_circular_cylinder::Set_max_step(int Step)
{
	_max_step = Step;
}

void flow_past_a_circular_cylinder::Set_lagrange_area()
{
	_lagrange_area_x1 = _center_x - _r - 4;

	_lagrange_area_x2 = _center_x + _r + 4;

	_lagrange_area_y1 = _center_y - _r - 4;

	_lagrange_area_y2 = _center_y + _r + 4;

	cout << "lagrange_area_x1=" << _lagrange_area_x1 << endl;

	cout << "lagrange_area_x2=" << _lagrange_area_x2 << endl;

	cout << "lagrange_area_y1=" << _lagrange_area_y1 << endl;

	cout << "lagrange_area_y2=" << _lagrange_area_y2 << endl;



}

void flow_past_a_circular_cylinder::Set_Solid()
{
	Solid = new Lagrange_field(_center_x,_center_y,_r,_Number_of_node);
}

void flow_past_a_circular_cylinder::Set_Fuild()
{
	Fluid = new euler_field(X_Max,Y_Max,_Re,_L,_u,_max_step);
}

void flow_past_a_circular_cylinder::Set_Interreaction()
{
	_interreaction = new Interreaction(X_Max, Y_Max, _Number_of_node);


#pragma omp parallel for
	for (int i = 0; i < X_Max + 1; i++)
	{
		for (int j = 0; j < Y_Max + 1; j++)
		{
			_interreaction->_Is_Compute[i][j] = false;

			double euler_x = i;
			double euler_y = j;


			for (int node = 0; node < _Number_of_node; node++)
			{
				double lagrange_x = Solid->lagrange_node[node]->Get_Position_x();
				double lagrange_y = Solid->lagrange_node[node]->Get_Position_y();

				double distance = Distance(euler_x, euler_y, lagrange_x, lagrange_y);
				if (distance <= 3) {
					_interreaction->_Is_Compute[i][j] = true;

				}
			}


			//cout << i << "," << j << "," << _interreaction->_Is_Compute[i][j] <<endl;

			/*

			if (_interreaction->_Is_Compute[i][j]) {
				//cout << i << "," << j << "," << _interreaction->_Is_Compute[i][j] << endl;
				cout << i << "\t"<< j<<";"<<endl;
			}

			*/
		
		}
	}







}

double flow_past_a_circular_cylinder::Get_X_Max()
{
	return X_Max;
}

double flow_past_a_circular_cylinder::Get_Y_Max()
{
	return Y_Max;
}

double flow_past_a_circular_cylinder::Get_center_x()
{
	return _center_x;
}

double flow_past_a_circular_cylinder::Get_center_y()
{
	return _center_y;
}

double flow_past_a_circular_cylinder::Get_r()
{
	return _r;
}

double flow_past_a_circular_cylinder::Get_Re()
{
	return _Re;
}

double flow_past_a_circular_cylinder::Get_L()
{
	return _L;
}

double flow_past_a_circular_cylinder::Get_number_of_threads()
{
	return _number_of_threads;
}

double flow_past_a_circular_cylinder::Get_u()
{
	return _u;
}

Lagrange_field * flow_past_a_circular_cylinder::Get_Solid()
{
	return Solid;
}

euler_field * flow_past_a_circular_cylinder::Get_Fuild()
{
	return Fluid;
}


//explicit diffuse direct_forcing IB_LBM
//step a
void flow_past_a_circular_cylinder::u_ij_step()
{
#pragma omp parallel for
	for (int i = 0; i < Fluid->euler_node.size(); i++)
	{
		for (int j = 0; j <  Fluid->euler_node[0].size(); j++)
		{
			double u = Fluid->euler_node[i][j]->Get_Velocity_x();
			double v = Fluid->euler_node[i][j]->Get_Velocity_y();

			Fluid->euler_node[i][j]->Set_u_v_noF(u, v);
		}
	}
}

//step b
void flow_past_a_circular_cylinder::Unforce_velocity_interpolation()
{
#pragma omp parallel for
	for (signed int node = 0; node < Solid->lagrange_node.size(); node++)
	{
		double ub_noF = 0;
		double vb_noF = 0;
		for (signed int i = _lagrange_area_x1; i < _lagrange_area_x2; i++)
		{
			for (signed int j = _lagrange_area_y1; j < _lagrange_area_y2; j++)
			{
				double h = 1;
				double D = D_function(Fluid->euler_node[i][j]->Get_Position_x(), Fluid->euler_node[i][j]->Get_Position_y(),
					Solid->lagrange_node[node]->Get_Position_x(), Solid->lagrange_node[node]->Get_Position_y());

				ub_noF += Fluid->euler_node[i][j]->Get_u_noF() * D * h * h;
				vb_noF += Fluid->euler_node[i][j]->Get_v_noF() * D * h * h;

			}
		}
		Solid->lagrange_node[node]->Set_ub_vb_noF(ub_noF, vb_noF);

	}
}

//step c
void flow_past_a_circular_cylinder::Boundary_force_evaluation_on_Xb()
{
#pragma omp parallel for
	for (signed int i = 0; i < Solid->lagrange_node.size(); i++)
	{
		double Fx = 0;
		double Fy = 0;
		//rho 为平均密度。。。这里没写完
		double rho = 1;
		double delta_t = 1;
		Fx = 2 * rho * (Solid->lagrange_node[i]->Get_Ub() - Solid->lagrange_node[i]->get_ub_nof()) / delta_t;
		Fy = 2 * rho * (Solid->lagrange_node[i]->Get_Vb() - Solid->lagrange_node[i]->get_vb_nof()) / delta_t;
		Solid->lagrange_node[i]->Set_F(Fx, Fy);
	}

}

//step d
void flow_past_a_circular_cylinder::Force_distribution_on_Xij()
{
#pragma omp parallel for
	for (signed int i = _lagrange_area_x1; i < _lagrange_area_x2; i++)
	{
		for (signed int j = _lagrange_area_y1; j < _lagrange_area_y2; j++)
		{

			double fx = 0;
			double fy = 0;
			for (signed int node = 0; node < Solid->lagrange_node.size(); node++)
			{
				double x1 = Solid->lagrange_node[node]->Get_Position_x();
				double y1 = Solid->lagrange_node[node]->Get_Position_y();

				double x2 = Fluid->euler_node[i][j]->Get_Position_x();
				double y2 = Fluid->euler_node[i][j]->Get_Position_y();

				double delta_sb;
				delta_sb = Solid->Get_delta_Sb();

				fx += Solid->lagrange_node[node]->Get_Fx() * D_function(x1, y1, x2, y2) * delta_sb;
				fy += Solid->lagrange_node[node]->Get_Fy() * D_function(x1, y1, x2, y2) * delta_sb;

			}
			Fluid->euler_node[i][j]->Set_Body_force(fx, fy);

		}
	}
}

//step e
void flow_past_a_circular_cylinder::Update_of_velocity_on_xij()
{
	double delta_t = 1;
	double rho = 1;
	#pragma omp parallel for
	for (signed int i = 0; i < Fluid->euler_node.size(); i++)
	{
		for (signed int j = 0; j < Fluid->euler_node[0].size(); j++)
		{
			double u = Fluid->euler_node[i][j]->Get_u_noF() + delta_t / 2 / rho * Fluid->euler_node[i][j]->Get_Body_force_fx();
			double v = Fluid->euler_node[i][j]->Get_v_noF() + delta_t / 2 / rho * Fluid->euler_node[i][j]->Get_Body_force_fy();

			Fluid->euler_node[i][j]->Set_Velocity(u, v);

		}
	}




}

void flow_past_a_circular_cylinder::First_forcing_step()
{
	u_ij_step();

	Unforce_velocity_interpolation();

	Boundary_force_evaluation_on_Xb();

	Force_distribution_on_Xij();


	Update_of_velocity_on_xij();
}

void flow_past_a_circular_cylinder::Collision_step()
{
	Fluid->Collision();

}

void flow_past_a_circular_cylinder::Secord_forcing_step()
{
#pragma omp parallel for
	for (signed int i = 0; i < Fluid->euler_node.size(); i++)
	{
		for (signed int j = 0; j < Fluid->euler_node[0].size(); j++)
		{


			for (signed int k = 0; k < 9; k++)
			{
				double F = force_distribution_function(Fluid->euler_node[i][j]->Get_Body_force_fx(), Fluid->euler_node[i][j]->Get_Body_force_fy(),
					Fluid->euler_node[i][j]->Get_Velocity_x(), Fluid->euler_node[i][j]->Get_Velocity_y(), k);
				
				Fluid->euler_node[i][j]->Set_F(F, k);


			}
		}
	}
#pragma omp parallel for
	for (signed int i = 0; i < Fluid->euler_node.size(); i++)
	{
		for (signed int j = 0; j < Fluid->euler_node[0].size(); j++)
		{

			for (signed int k = 0; k < 9; k++)
			{
				double fk = Fluid->euler_node[i][j]->Get_f(k) + Fluid->euler_node[i][j]->Get_F(k);
				Fluid->euler_node[i][j]->Set_f(fk, k);
			}
		}
	}


}

void flow_past_a_circular_cylinder::Stream_step()
{
	Fluid->Stream();
}

void flow_past_a_circular_cylinder::Macroscopic()
{
	Fluid->Macroscopic();
	Fluid->Boundary_Macroscopic();
}

void flow_past_a_circular_cylinder::evolution()
{
	Fluid->Boundary_Condition();
	Fluid->Output_grid();
	Output_Parameter();
	for (int i = 0; i < 10001; i++)
	{
		First_forcing_step();
		Collision_step();
		Secord_forcing_step();
		Stream_step();
		Macroscopic();

		if (i % 2500 == 0) {
			Output_Fluid(_Re, i);
			Output_Solid(_Re, i);
			cout << "evolution:"<<i << endl;
		}

	}
}

void flow_past_a_circular_cylinder::Set_up_interreaction_data()
{
	Set_Interreaction();

	//euler node sent massage
	_interreaction->Set_euler_node_number();

	_interreaction->init_euler_data();

	euler_node_to_interreaction();


	//lagrange node sent massage

	_interreaction->Set_lagrange_node_number();

	_interreaction->init_lagrange_data();

	lagrange_node_to_interreation();


	_interreaction->Set_valid_cal();

	_interreaction->Evolution_Step();

	interreaction_to_euler_node();
	interreaction_to_lagrange_node();


}

void flow_past_a_circular_cylinder::multi_evolution()
{

	Fluid->Boundary_Condition();
	Fluid->Output_grid();
	Output_Parameter();
	//Set_up_interreaction_data();
	//_interreaction->Output_euler_data();
	for (int i = 0; i < _max_step; i++)
	{
		multi_direct_forcing_step();

		Collision_step();
		Secord_forcing_step();
		Stream_step();
		Macroscopic();

		if (i % 5000 == 0) {
			cout << i << endl;
			Output_Fluid(_Re, i+1);
			Output_Solid(_Re, i+1);
			//_interreaction->output_euler_node(_Re, i+1);
		}


		//cout << i << endl;
		//Output_Fluid(_Re, _max_step + i);
		//Output_Solid(_Re, _max_step + i);

		
		
		

		//delete _interreaction;

	}
	

	//_interreaction->Output_euler_data();
	//_interreaction->Output_laghrange_data();

	//multi_direct_forcing_step();


	//Output_Fluid(_Re, 999);
	//Output_Solid(_Re, 999);


}

void flow_past_a_circular_cylinder::multi_direct_forcing_step()
{

	Set_Interreaction();

	//euler node sent massage
	_interreaction->Set_euler_node_number();

	_interreaction->init_euler_data();

	euler_node_to_interreaction();
	

	//lagrange node sent massage
	
	_interreaction->Set_lagrange_node_number();

	_interreaction->init_lagrange_data();
	
	lagrange_node_to_interreation();

	


	_interreaction->Set_valid_cal();

	_interreaction->Evolution_Step();
	
	interreaction_to_euler_node();
	interreaction_to_lagrange_node();
	
	//cout <<"---------------------------"<<endl;
	//interreaction_to_lagrange_node();
	//interreaction_to_euler_node();


	//output

	//_interreaction->Output_euler_data();
	//_interreaction->Output_laghrange_data();


	delete _interreaction;
	

}

void flow_past_a_circular_cylinder::euler_node_to_interreaction()
{
#pragma omp parallel for
	for (int i = 0; i < _interreaction->_Euler_data.size(); i++)
	{
		int x = _interreaction->_Euler_data[i]->Get_number_x();
		int y = _interreaction->_Euler_data[i]->Get_number_y();


		//euler_data
		double u = Fluid->euler_node[x][y]->Get_Velocity_x();
		double v = Fluid->euler_node[x][y]->Get_Velocity_y();

		double x_position = Fluid->euler_node[x][y]->Get_Position_x();
		double y_position = Fluid->euler_node[x][y]->Get_Position_y();


		//interreaction_data

		_interreaction->_Euler_data[i]->Set_u_ij_m(u, v, 0);
		_interreaction->_Euler_data[i]->Set_Position(x_position, y_position);

		//_interreaction->Evolution_Step();


	}


}

void flow_past_a_circular_cylinder::lagrange_node_to_interreation()
{
#pragma omp parallel for
	for (int i = 0; i < Solid->lagrange_node.size(); i++)
	{
		double x = Solid->lagrange_node[i]->Get_Position_x();
		double y = Solid->lagrange_node[i]->Get_Position_y();

		double U_b_x = Solid->lagrange_node[i]->Get_Ub();
		double U_b_y = Solid->lagrange_node[i]->Get_Vb();


		double delta_sb = Solid->lagrange_node[i]->Get_delta_sb();




		//cout << delta_sb;




		_interreaction->_lagrange_data[i]->Set_grid_xy(x,y);

		_interreaction->_lagrange_data[i]->Set_delta_sb(delta_sb);

		_interreaction->_lagrange_data[i]->Set_U_b(U_b_x, U_b_y);
		//_interreaction->_lagrange_data[i]->Set_U_b(1,2 );


	}



}

void flow_past_a_circular_cylinder::interreaction_to_lagrange_node()
{

#pragma omp parallel for
	for (int i = 0; i < _interreaction->_lagrange_data.size(); i++)
	{
		double node_number = _interreaction->_lagrange_data[i]->Get_node_number();

		//cout << node_number << endl;

		double fbx = _interreaction->_lagrange_data[i]->Get_result_F_b_x();
		double fby = _interreaction->_lagrange_data[i]->Get_result_F_b_y();


		//cout << fx << "\t" <<fy << endl;

		Solid->lagrange_node[i]->Set_F(fbx, fby);



	}



}

void flow_past_a_circular_cylinder::interreaction_to_euler_node()
{
#pragma omp parallel for
	for (int i = 0; i < _interreaction->_Euler_data.size(); i++)
	{
		double x = _interreaction->_Euler_data[i]->Get_number_x();
		double y = _interreaction->_Euler_data[i]->Get_number_y();

		double ux = _interreaction->_Euler_data[i]->Get_result_u_ij_x();
		double uy = _interreaction->_Euler_data[i]->Get_result_u_ij_y();

		double fx = _interreaction->_Euler_data[i]->Get_result_F_ij_x();
		double fy = _interreaction->_Euler_data[i]->Get_result_F_ij_y();


		
		//cout << x << "\t" << y << endl;
		//cout << ux << "\t" << uy << endl;
		//cout << fx << "\t" << fy << endl;





		Fluid->euler_node[x][y]->Set_Body_force(fx, fy);
		Fluid->euler_node[x][y]->Set_Velocity(ux, uy);

/*
		cout << x << "\t" << y << endl;

		cout << Fluid->euler_node[x][y]->Get_Velocity_x() << "\t" << Fluid->euler_node[x][y]->Get_Velocity_y();
		cout << endl;

		cout << Fluid->euler_node[x][y]->Get_Body_force_fx() << "\t" << Fluid->euler_node[x][y]->Get_Body_force_fy();


		cout << endl;
*/



	}







}




void flow_past_a_circular_cylinder::Output_Parameter()
{
	cout << "Grid_X = " << X_Max << endl;
	cout << "Grid_Y = " << Y_Max << endl;



	cout <<"number_of_node="<<_Number_of_node<< endl;
	cout << "center = " << "(" <<_center_x<<","<< _center_y<< ")" << endl;
	cout << "r = " << _r<<endl;
	cout << endl;

	cout << "Re = " << _Re << endl;
	cout << "u = " << _u << endl;
	cout << "L = " << _L << endl;
	cout <<"number_of_threads = " <<omp_get_thread_num()<<endl;
	cout << "max_step = " << _max_step<<endl;

}

void flow_past_a_circular_cylinder::Output_Data()
{
}

void flow_past_a_circular_cylinder::Output_Solid(double value, int step)
{
	Solid->Output_Solid(value, step);
	
}

void flow_past_a_circular_cylinder::Output_Fluid(double value, int step)
{
	Fluid->Output_field(value, step); 
	Fluid->output_all(value, step);

}

flow_past_a_circular_cylinder::~flow_past_a_circular_cylinder()
{





}

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



//mluti method
flow_past_a_circular_cylinder::flow_past_a_circular_cylinder(int Grid_X, int Grid_Y,
	double center_x, double center_y, double r, int number_of_node, 
	double Re, double L, double u, int number_of_threads, 
	int max_step, int NF)
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

	Set_NF(NF);
	//cout <<"nf = "<< _NF;
	Set_Solid(_NF);
	Set_Fuild(_NF);



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

	//cout << "lagrange_area_x1=" << _lagrange_area_x1 << endl;

	//cout << "lagrange_area_x2=" << _lagrange_area_x2 << endl;

	//cout << "lagrange_area_y1=" << _lagrange_area_y1 << endl;

	//cout << "lagrange_area_y2=" << _lagrange_area_y2 << endl;



}


void flow_past_a_circular_cylinder::Set_NF(int nf_number)
{
	_NF = nf_number;
}

void flow_past_a_circular_cylinder::Set_Solid()
{
	Solid = new Lagrange_field(_center_x,_center_y,_r,_Number_of_node);
}

void flow_past_a_circular_cylinder::Set_Fuild()
{
	Fluid = new euler_field(X_Max,Y_Max,_Re,_L,_u,_max_step);
}

void flow_past_a_circular_cylinder::Set_Solid(int nf)
{
	Solid = new Lagrange_field(_center_x, _center_y, _r, _Number_of_node,_NF);
}

void flow_past_a_circular_cylinder::Set_Fuild(int nf)
{
	Fluid = new euler_field(X_Max, Y_Max, _Re, _L, _u, _max_step, _NF);
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

int flow_past_a_circular_cylinder::Get_NF()
{
	return _NF;
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

		int delta_function_area = 3;
		int x_min = int(Solid->lagrange_node[node]->Get_Position_x()) - delta_function_area;
		int x_max = int(Solid->lagrange_node[node]->Get_Position_x()) + delta_function_area;
		int y_min = int(Solid->lagrange_node[node]->Get_Position_y()) - delta_function_area;
		int y_max = int(Solid->lagrange_node[node]->Get_Position_y()) + delta_function_area;

		for (signed int i = x_min; i < x_max; i++)
		{
			for (signed int j = y_min; j < y_max; j++)
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

	/*
	#pragma omp parallel for
	for (signed int i = 0; i < Fluid->euler_node.size(); i++)
	{
		for (signed int j = 0; j < Fluid->euler_node[0].size(); j++)
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
	*/
	
	//加速版


	//清空Fij

	
#pragma omp parallel for
	for (signed int i = _lagrange_area_x1; i < _lagrange_area_x2; i++)
	{
		for (signed int j = _lagrange_area_y1; j < _lagrange_area_y2; j++)
		{

			double fx = 0;
			double fy = 0;
			
			Fluid->euler_node[i][j]->Set_Body_force(fx, fy);

		}
	}
	



	//赋值累加力项
	#pragma omp parallel for
	for (signed int node = 0; node < Solid->lagrange_node.size(); node++)
	{
		

		int delta_function_area = 3;
		int x_min = int(Solid->lagrange_node[node]->Get_Position_x()) - delta_function_area;
		int x_max = int(Solid->lagrange_node[node]->Get_Position_x()) + delta_function_area;
		int y_min = int(Solid->lagrange_node[node]->Get_Position_y()) - delta_function_area;
		int y_max = int(Solid->lagrange_node[node]->Get_Position_y()) + delta_function_area;

		for (signed int i = x_min; i < x_max; i++)
		{
			for (signed int j = y_min; j < y_max; j++)
			{
				double fx = Fluid->euler_node[i][j]->Get_Body_force_fx();
				double fy = Fluid->euler_node[i][j]->Get_Body_force_fy();


				double x1 = Solid->lagrange_node[node]->Get_Position_x();
				double y1 = Solid->lagrange_node[node]->Get_Position_y();

				double x2 = Fluid->euler_node[i][j]->Get_Position_x();
				double y2 = Fluid->euler_node[i][j]->Get_Position_y();

				double delta_sb;
				delta_sb = Solid->Get_delta_Sb();

				fx += Solid->lagrange_node[node]->Get_Fx() * D_function(x1, y1, x2, y2) * delta_sb;
				fy += Solid->lagrange_node[node]->Get_Fy() * D_function(x1, y1, x2, y2) * delta_sb;

				Fluid->euler_node[i][j]->Set_Body_force(fx, fy);

				
			}
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

void flow_past_a_circular_cylinder::mluti_method_step_a()
{




#pragma omp parallel for
	for (int i = 0; i < Fluid->euler_node.size(); i++)
	{
		for (int j = 0; j < Fluid->euler_node[0].size(); j++)
		{
			double u = Fluid->euler_node[i][j]->Get_Velocity_x();
			double v = Fluid->euler_node[i][j]->Get_Velocity_y();

			Fluid->euler_node[i][j]->Set_NF_u(0, u, v);

			//Fluid->euler_node[i][j]->Set_Velocity(u, v);
		}
	}





}
//
//void flow_past_a_circular_cylinder::mluti_method_step_b()
//{
//#pragma omp parallel for
//	for (signed int node = 0; node < Solid->lagrange_node.size(); node++)
//	{
//		double ub = 0;
//		double vb = 0;
//
//		int delta_function_area = 3;
//		int x_min = int(Solid->lagrange_node[node]->Get_Position_x()) - delta_function_area;
//		int x_max = int(Solid->lagrange_node[node]->Get_Position_x()) + delta_function_area;
//		int y_min = int(Solid->lagrange_node[node]->Get_Position_y()) - delta_function_area;
//		int y_max = int(Solid->lagrange_node[node]->Get_Position_y()) + delta_function_area;
//
//		for (signed int i = x_min; i < x_max; i++)
//		{
//			for (signed int j = y_min; j < y_max; j++)
//			{
//				double h = 1;
//				double D = D_function(Fluid->euler_node[i][j]->Get_Position_x(), Fluid->euler_node[i][j]->Get_Position_y(),
//					Solid->lagrange_node[node]->Get_Position_x(), Solid->lagrange_node[node]->Get_Position_y());
//
//				ub += Fluid->euler_node[i][j]->Get_NF_u(0) * D * h * h;
//				vb += Fluid->euler_node[i][j]->Get_NF_v(0) * D * h * h;
//
//			}
//		}
//		Solid->lagrange_node[node]->Set_NF_ub_vb(ub, vb, 0);
//
//	}
//
//
//}




void flow_past_a_circular_cylinder::mluti_method_step_b(){
	
#pragma omp parallel for
	for (signed int node = 0; node < Solid->lagrange_node.size(); node++)
	{
		double ub_noF = 0;
		double vb_noF = 0;

		int delta_function_area = 3;
		int x_min = int(Solid->lagrange_node[node]->Get_Position_x()) - delta_function_area;
		int x_max = int(Solid->lagrange_node[node]->Get_Position_x()) + delta_function_area;
		int y_min = int(Solid->lagrange_node[node]->Get_Position_y()) - delta_function_area;
		int y_max = int(Solid->lagrange_node[node]->Get_Position_y()) + delta_function_area;

		for (signed int i = x_min; i < x_max; i++)
		{
			for (signed int j = y_min; j < y_max; j++)
			{
				double h = 1;
				double D = D_function(Fluid->euler_node[i][j]->Get_Position_x(), Fluid->euler_node[i][j]->Get_Position_y(),
					Solid->lagrange_node[node]->Get_Position_x(), Solid->lagrange_node[node]->Get_Position_y());

				ub_noF += Fluid->euler_node[i][j]->Get_u_noF() * D * h * h;
				vb_noF += Fluid->euler_node[i][j]->Get_v_noF() * D * h * h;

			}
		}

		
		Solid->lagrange_node[node]->Set_ub_vb_noF(ub_noF, vb_noF);



		Solid->lagrange_node[node]->Set_NF_ub_vb(ub_noF, vb_noF,0);

	}



}


void flow_past_a_circular_cylinder::mluti_method_step_c(int _m)
{
#pragma omp parallel for
	for (signed int i = 0; i < Solid->lagrange_node.size(); i++)
	{
		double Fx_m = 0;
		double Fy_m = 0;
		

		double rho = 0.7 ;
		double delta_t = 1;
		//Fx = 2 * rho * (Solid->lagrange_node[i]->Get_Ub() - Solid->lagrange_node[i]->get_ub_nof()) / delta_t;
		//Fy = 2 * rho * (Solid->lagrange_node[i]->Get_Vb() - Solid->lagrange_node[i]->get_vb_nof()) / delta_t;
		//Solid->lagrange_node[i]->Set_F(Fx, Fy);
		
		double Ub = Solid->lagrange_node[i]->Get_Ub();
		double Vb = Solid->lagrange_node[i]->Get_Vb();

		double ub = Solid->lagrange_node[i]->Get_NF_ub(_m - 1);
		double vb = Solid->lagrange_node[i]->Get_NF_vb(_m - 1);

		Fx_m = 2 * rho * (Ub - ub) / delta_t;
		Fy_m = 2 * rho * (Vb - vb) / delta_t;
		
		Solid->lagrange_node[i]->Set_NF_F(Fx_m, Fy_m, _m);
	
	}




}

void flow_past_a_circular_cylinder::mluti_method_step_d(int _m)
{
	
	//clear nf Fij
//#pragma omp parallel for
	/*for (signed int i = _lagrange_area_x1; i < _lagrange_area_x2; i++)
	{
		for (signed int j = _lagrange_area_y1; j < _lagrange_area_y2; j++)
		{

			double fx = 0;
			double fy = 0;

			Fluid->euler_node[i][j]->Set_NF_F(_m, fx, fy);

		}
	}
	*/


	/*
#pragma omp parallel for
	for (signed int node = 0; node < Solid->lagrange_node.size(); node++)
	{


		int delta_function_area = 3;
		int x_min = int(Solid->lagrange_node[node]->Get_Position_x()) - delta_function_area;
		int x_max = int(Solid->lagrange_node[node]->Get_Position_x()) + delta_function_area;
		int y_min = int(Solid->lagrange_node[node]->Get_Position_y()) - delta_function_area;
		int y_max = int(Solid->lagrange_node[node]->Get_Position_y()) + delta_function_area;

		for (signed int i = x_min; i < x_max; i++)
		{
			for (signed int j = y_min; j < y_max; j++)
			{
				//double fx = Fluid->euler_node[i][j]->Get_Body_force_fx();
				//double fy = Fluid->euler_node[i][j]->Get_Body_force_fy();

				double fx = Fluid->euler_node[i][j]->Get_NF_Fx(_m);
				double fy = Fluid->euler_node[i][j]->Get_NF_Fy(_m);


				double x1 = Solid->lagrange_node[node]->Get_Position_x();
				double y1 = Solid->lagrange_node[node]->Get_Position_y();

				double x2 = Fluid->euler_node[i][j]->Get_Position_x();
				double y2 = Fluid->euler_node[i][j]->Get_Position_y();

				double delta_sb;
				delta_sb = Solid->Get_delta_Sb();




				fx += Solid->lagrange_node[node]->Get_NF_Fx(_m) * D_function(x1, y1, x2, y2) * delta_sb;
				fy += Solid->lagrange_node[node]->Get_NF_Fy(_m) * D_function(x1, y1, x2, y2) * delta_sb;

				Fluid->euler_node[i][j]->Set_NF_F(_m, fx, fy);


			}

		}


	}

	*/

	#pragma omp parallel for
	for (int i = 0; i < Fluid->euler_node.size(); i++)
	{
		for (int j = 0; j < Fluid->euler_node[0].size(); j++)
		{
			double Fx = 0;
			double Fy = 0;

			for (int node = 0; node < Solid->lagrange_node.size(); node++)
			{
				double x1 = Solid->lagrange_node[node]->Get_Position_x();
				double y1 = Solid->lagrange_node[node]->Get_Position_y();

				double x2 = Fluid->euler_node[i][j]->Get_Position_x();
				double y2 = Fluid->euler_node[i][j]->Get_Position_y();

				if (D_function(x1, y1, x2, y2) != 0)
				{

\
					
					double delta_sb;
					delta_sb = Solid->Get_delta_Sb();


					/*cout << "x = " << Fluid->euler_node[i][j]->Get_Position_x() << "\t";
					cout << "y = " << Fluid->euler_node[i][j]->Get_Position_y() << "\t";
					cout << "fx= " << Solid->lagrange_node[node]->Get_NF_Fx(_m)* D_function(x1, y1, x2, y2) * delta_sb << "\t";
					cout << "fy= " << Solid->lagrange_node[node]->Get_NF_Fy(_m)* D_function(x1, y1, x2, y2) * delta_sb << "\t";*/


					/*cout << endl;*/
					

					Fx += Solid->lagrange_node[node]->Get_NF_Fx(_m)* D_function(x1, y1, x2, y2) * delta_sb;
					Fy += Solid->lagrange_node[node]->Get_NF_Fy(_m)* D_function(x1, y1, x2, y2) * delta_sb;
					
				}
	
			}

			////输出部分
			//if (_step == 300)
			//{
			//	
			//	if (Fx != 0 || Fy != 0)
			//	{
			//		cout << "before set" << endl;
			//		cout << "stepe_1" << endl;
			//		cout << " m =" << _m << endl;
			//		cout << "x = " << Fluid->euler_node[i][j]->Get_Position_x() << "\t";
			//		cout << "y = " << Fluid->euler_node[i][j]->Get_Position_y() << "\t";

			//		cout <<Fx<< "\t";
			//		cout <<Fy<< "\t";
			//		cout << endl;
			//	}
			//}

			Fluid->euler_node[i][j]->Set_NF_F(_m, Fx, Fy);

			//if (Fx !=0 || Fy !=0)
			//{
			//	cout << "x = " << i << "\t";
			//	cout << "y = " << j << "\t";
			//	cout << "Fx = " << Fx << "\t";
			//	cout << "Fy = " << Fy << "\t";


			//	cout << Fluid->euler_node[i][j]->Get_NF_Fx(_m) << "\t";

			//	cout << Fluid->euler_node[i][j]->Get_NF_Fy(_m) << "\t";


			//	cout << endl;

			//}

			

			
		/*	if (_step == 300)
			{
				
				if (Fx != 0 || Fy != 0)
				{
					cout << "after set" << endl;
					cout << "stepe_2" << endl;
					cout << " m =" << _m << endl;

					cout << "x = " << Fluid->euler_node[i][j]->Get_Position_x() << "\t";
					cout << "y = " << Fluid->euler_node[i][j]->Get_Position_y() << "\t";

					cout << Fluid->euler_node[i][j]->Get_NF_Fx(_m) << "\t";
					cout << Fluid->euler_node[i][j]->Get_NF_Fy(_m) << "\t";
					cout << endl;

				}
			}*/




		}
	}







}

void flow_past_a_circular_cylinder::mluti_method_step_e(int _m)
{
	double delta_t = 1;
	double rho = 0.5;
	#pragma omp parallel for
	for (signed int i = 0; i < Fluid->euler_node.size(); i++)
	{
		for (signed int j = 0; j < Fluid->euler_node[0].size(); j++)
		{
			//double u = Fluid->euler_node[i][j]->Get_u_noF() + delta_t / 2 / rho * Fluid->euler_node[i][j]->Get_Body_force_fx();
			//double v = Fluid->euler_node[i][j]->Get_v_noF() + delta_t / 2 / rho * Fluid->euler_node[i][j]->Get_Body_force_fy();
			
			double u = Fluid->euler_node[i][j]->Get_NF_u(_m - 1) + delta_t / 2 / rho * Fluid->euler_node[i][j]->Get_NF_Fx(_m);
			double v = Fluid->euler_node[i][j]->Get_NF_v(_m - 1) + delta_t / 2 / rho * Fluid->euler_node[i][j]->Get_NF_Fy(_m);



			//Fluid->euler_node[i][j]->Set_Velocity(u, v);

			Fluid->euler_node[i][j]->Set_NF_u(_m, u, v);

		}
	}




}

void flow_past_a_circular_cylinder::mluti_method_step_f(int _m)
{
#pragma omp parallel for
	for (signed int node = 0; node < Solid->lagrange_node.size(); node++)
	{
		double ub = 0;
		double vb = 0;

		int delta_function_area = 3;
		int x_min = int(Solid->lagrange_node[node]->Get_Position_x()) - delta_function_area;
		int x_max = int(Solid->lagrange_node[node]->Get_Position_x()) + delta_function_area;
		int y_min = int(Solid->lagrange_node[node]->Get_Position_y()) - delta_function_area;
		int y_max = int(Solid->lagrange_node[node]->Get_Position_y()) + delta_function_area;

		for (signed int i = x_min; i < x_max; i++)
		{
			for (signed int j = y_min; j < y_max; j++)
			{
				double h = 1;
				double D = D_function(Fluid->euler_node[i][j]->Get_Position_x(), Fluid->euler_node[i][j]->Get_Position_y(),
					Solid->lagrange_node[node]->Get_Position_x(), Solid->lagrange_node[node]->Get_Position_y());

				ub += Fluid->euler_node[i][j]->Get_NF_u(_m) * D * h * h;
				vb += Fluid->euler_node[i][j]->Get_NF_v(_m) * D * h * h;

			}
		}
		Solid->lagrange_node[node]->Set_NF_ub_vb(ub, vb, _m);

	}


}

void flow_past_a_circular_cylinder::mluti_method_step_g()
{
	//Fb
	for (int i = 0; i < Solid->lagrange_node.size(); i++)
	{
		double Fx = 0;
		double Fy = 0;
		for (int m = 0; m < _NF; m++)
		{

			Fx += Solid->lagrange_node[i]->Get_NF_Fx(m);
			Fy += Solid->lagrange_node[i]->Get_NF_Fy(m);

		}


		Solid->lagrange_node[i]->Set_F(Fx, Fy);

	}




	//ub

	/*for (int i = 0; i < Solid->lagrange_node.size(); i++)
	{
		double ub = Solid->lagrange_node[i]->Get_NF_ub(_NF-1);
		double vb = Solid->lagrange_node[i]->Get_NF_vb(_NF-1);

		Solid->lagrange_node[i]->Set_ub_vb_noF(ub,vb);
	}


*/


	//F_ij

	for (int i = 0; i < Fluid->euler_node.size(); i++)
	{
		for (int j = 0; j < Fluid->euler_node[0].size(); j++)
		{
			double Fx = 0;
			double Fy = 0;

			for (int m = 0; m < _NF; m++)
			{
				Fx += Fluid->euler_node[i][j]->Get_NF_Fx(m);
				Fy += Fluid->euler_node[i][j]->Get_NF_Fy(m);
			}

			Fluid->euler_node[i][j]->Set_Body_force(Fx, Fy);


		}
	}






	//uij

	for (int i = 0; i < Fluid->euler_node.size(); i++)
	{
		for (int j = 0; j < Fluid->euler_node[0].size(); j++)
		{
			double u = Fluid->euler_node[i][j]->Get_NF_u(_NF - 1);
			double v = Fluid->euler_node[i][j]->Get_NF_v(_NF - 1);

			Fluid->euler_node[i][j]->Set_Velocity(u, v);

		}
	}








}






void flow_past_a_circular_cylinder::mluti_method_first_forcing_step()
{

	//clear_NF_Parameter();
	mluti_method_step_a();	
	mluti_method_step_b();


	for (int m = 1; m < _NF; m++)
	{
		mluti_method_step_c(m);

		mluti_method_step_d(m);

		mluti_method_step_e(m);

		mluti_method_step_f(m);

	}
	mluti_method_step_g();


}

void flow_past_a_circular_cylinder::mluti_method_evolution()
{
	
	Fluid->Boundary_Condition();
	Fluid->Output_grid();
	Output_Parameter();
	for (int i = 0; i < _max_step; i++)
	{

		mluti_method_first_forcing_step();
		Collision_step();
		Secord_forcing_step();
		Stream_step();
		Macroscopic();

		/*	Solid->Output_Solid_NF(0, 0, i);
			Fluid->Output_Fluid_NF(0, 0, i);
			Output_Fluid(_Re, i);
			Output_Solid(_Re, i);*/

		if (i % 1000 == 0) {
			Output_Fluid(_Re, i);
			Output_Solid(_Re, i);

			Fluid->Output_all(_Re, i);
			Solid->Output_Solid_NF(0, 0, i);
			Fluid->Output_Fluid_NF(0, 0, i);

			//Output_mluti_test();
	


			cout << i << endl;
		}
		_step++;

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
	for (int i = 0; i < _max_step; i++)
	{
		First_forcing_step();
		Collision_step();
		Secord_forcing_step();
		Stream_step();
		Macroscopic();

		if (i % 1000 == 0) {
			Output_Fluid(_Re, i);
			Fluid->Output_all(_Re, i);

			Output_Solid(_Re, i);
			cout << i << endl;
		}

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
	//cout << "max_step = " << Solid->

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

}



void flow_past_a_circular_cylinder::mluti_test_step_a()
{
	//int m = 0;

	
	mluti_method_step_a();
	//Fluid->Output_Fluid_NF(0, 0, 20000);

}

void flow_past_a_circular_cylinder::mluti_test_step_b()
{
	mluti_method_step_b();
	//Solid->Output_Solid_NF(0, 0, 20000);

}

void flow_past_a_circular_cylinder::mluti_test_step_c(int m)
{

	mluti_method_step_c(m);
	
	//Fluid->Output_Fluid_NF(0, 0, 20000);
	//Solid->Output_Solid_NF(0, 0, 20000);

}

void flow_past_a_circular_cylinder::mluti_test_step_d(int m)
{
	mluti_method_step_d(m);

	//Fluid->Output_Fluid_NF(0, 0, 20000);
	//Solid->Output_Solid_NF(0, 0, 20000);


	//cout << "..........................................." << endl;
	//for (int i = 0; i < Fluid->euler_node.size(); i++)
	//{
	//	for (int j = 0; j < Fluid->euler_node.size(); j++)
	//	{
	//		
	//		double Fx = Fluid->euler_node[i][j]->Get_NF_Fx(m);
	//		double Fy = Fluid->euler_node[i][j]->Get_NF_Fy(m);

	//		if (Fx != 0 || Fy != 0)
	//		{
	//			cout << "x = " << i << "\t";
	//			cout << "y = " << j << "\t";
	//			cout << "Fx = " << Fx << "\t";
	//			cout << "Fy = " << Fy << "\t";


	//			cout << Fluid->euler_node[i][j]->Get_NF_Fx(m) << "\t";

	//			cout << Fluid->euler_node[i][j]->Get_NF_Fy(m) << "\t";


	//			cout << endl;

	//		}


	//	}
	//}



}

void flow_past_a_circular_cylinder::mluti_test_step_e(int m)
{
	
	mluti_method_step_e(m);
	//Fluid->Output_Fluid_NF(0, 0, 20000);
	//Solid->Output_Solid_NF(0, 0, 20000);

}

void flow_past_a_circular_cylinder::mluti_test_step_f(int m)
{
	mluti_method_step_f(m);
	//Fluid->Output_Fluid_NF(0, 0, 20000);
	//Solid->Output_Solid_NF(0, 0, 20000);


}

void flow_past_a_circular_cylinder::mluti_test_step_g()
{


	//Fluid->Output_Fluid_NF(0, 0, 20000);
	//Solid->Output_Solid_NF(0, 0, 20000);


	//Output_Fb();



	mluti_method_step_g();

	//Output_Fb();


}

void flow_past_a_circular_cylinder::clear_NF_Parameter()
{


	//Solid->Clear_NF_Fb();
	//Solid->Clear_NF_ub();

	//Fluid->Clear_NF_F_ij();
	//Fluid->Clear_NF_u();
}

void flow_past_a_circular_cylinder::Output_mluti_test()
{
	Fluid->Output_Fluid_NF(0, 0, 20000);
	Solid->Output_Solid_NF(0, 0, 20000);
}

void flow_past_a_circular_cylinder::Output_Fb()
{




	for (int i = 0; i < Solid->lagrange_node.size(); i++)
	{
		cout << endl;

		cout << "x = " << Solid->lagrange_node[i]->Get_Position_x() << "\t";
		cout << "y = " << Solid->lagrange_node[i]->Get_Position_y() << "\t";

		cout << "getfx = " << Solid->lagrange_node[i]->Get_Fx() << "\t";
		cout << "getfy = " << Solid->lagrange_node[i]->Get_Fy() << "\t";
		cout << endl;

	}


}









flow_past_a_circular_cylinder::~flow_past_a_circular_cylinder()
{





}

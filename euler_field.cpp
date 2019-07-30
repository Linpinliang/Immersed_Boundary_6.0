#include "euler_field.h"

using namespace std;

void euler_field::Set_Re(double Re)
{
	_Re = Re;
}

void euler_field::Set_L(double L)
{
	_L = L;
}

void euler_field::Set_tau()
{
	_tau = 0.5 + 3 * (_L * _u_inlet / _Re);

}

void euler_field::Set_u_inlet(double u)
{
	_u_inlet = u;
}

void euler_field::Set_Max_step(int N)
{
	_Max_Step = N;
	
}

void euler_field::Set_NF(int nf)
{
	_NF = nf;
}

double euler_field::Get_Re()
{
	return _Re;
}

double euler_field::Get_L()
{
	return _L;
}

double euler_field::Get_tau()
{
	return _tau;
}

double euler_field::Get_u_inlet()
{
	return _u_inlet;
}

int euler_field::Get_Max_step()
{
	return _Max_Step;
}

int euler_field::Get_NF()
{
	return 0;
}



euler_field::euler_field()
{
}

euler_field::euler_field(int x, int y, double Re, double L, double u_inlet,int _max_step)
{

	Nx = x;
	Ny = y;

	Set_Re(Re);
	Set_L(L);
	Set_u_inlet(u_inlet);
	Set_tau();
	Set_Max_step(_max_step);
	
	

	euler_node = vector<vector<Euler_Point*> >(Nx + 1, vector<Euler_Point*>(Ny + 1));

	for (int i = 0; i < Nx + 1; i++)
	{
		for (int j = 0; j < Ny + 1; j++)
		{
			Euler_Point* p1 = new Euler_Point(i, j);
			euler_node[i][j] = p1;
			//cout <<i<<"\t"<<j<<endl;
		}
	}
}

euler_field::euler_field(int x, int y, double Re, double L, double u_inlet, int _max_step, int nf_number)
{
	Nx = x;
	Ny = y;

	Set_Re(Re);
	Set_L(L);
	Set_u_inlet(u_inlet);
	Set_tau();
	Set_Max_step(_max_step);
	Set_NF(nf_number);
	cout << "eulerfiled nf =" << _NF<<endl;

	euler_node = vector<vector<Euler_Point*> >(Nx + 1, vector<Euler_Point*>(Ny + 1));

	for (int i = 0; i < Nx + 1; i++)
	{
		for (int j = 0; j < Ny + 1; j++)
		{
			Euler_Point* p1 = new Euler_Point(i, j,_NF);
			euler_node[i][j] = p1;
			//cout <<i<<"\t"<<j<<endl;
		}
	}

}





euler_field::~euler_field()
{
}



void euler_field::evolution()
{


	
	for (int i = 0; i < _Max_Step; i++)
	{
		Collision();
		Stream();
		Macroscopic();
		Boundary_Macroscopic();

		if (i%500==0)
		{
			Output_field(_Re, i);
			Output_grid();
		
		}


	}
	
}

void euler_field::Collision()
{
	#pragma omp parallel for
	for (signed int x = 1; x < Nx; x++)
	{
		for (signed int y = 1; y < Ny; y++)
		{
			for (int k = 0; k < Q; k++)
			{
				double f_k = euler_node[x][y]->Get_f(k);
				double f_k_eq = FeqStandardD2Q9(k, euler_node[x][y]->Get_rho(), euler_node[x][y]->Get_Velocity_x(), euler_node[x][y]->Get_Velocity_y());
				
				f_k = f_k - (f_k - f_k_eq) / _tau;
				euler_node[x][y]->Set_f(f_k, k);
			
			}

		}
	}

}

void euler_field::Stream()
{
	#pragma omp parallel for
	for (signed int x = 1; x < Nx; x++)
	{
		for (signed int y = 1; y < Ny; y++)
		{
			for (int k = 0; k < Q; k++)
			{
				int xp = x - e[k][0];
				int yp = y - e[k][1];

				euler_node[x][y]->Set_New_f(euler_node[xp][yp]->Get_f(k), k);

			}

		}
	}

}

void euler_field::Macroscopic()
{
	#pragma omp parallel for
	for (int x = 1; x < Nx; x++)
	{
		for (signed int y = 1; y < Ny; y++)
		{
			double u0 = euler_node[x][y]->Get_Velocity_x();
			double v0 = euler_node[x][y]->Get_Velocity_y();
			euler_node[x][y]->Set_Velocity0(u0, v0);

			euler_node[x][y]->Set_rho(0);
			euler_node[x][y]->Set_Velocity(0,0);


			double rho = 0;
			double u   = 0;
			double v   = 0;


			for (int k = 0; k < Q; k++)
			{
				double f = euler_node[x][y]->Get_New_f(k);
				euler_node[x][y]->Set_f(f,k);

				rho += euler_node[x][y]->Get_f(k);
				u += e[k][0] * euler_node[x][y]->Get_f(k);
				v += e[k][1] * euler_node[x][y]->Get_f(k);
				
			}

			u /= rho;
			v /= rho;

			euler_node[x][y]->Set_rho(rho);
			euler_node[x][y]->Set_Velocity(u, v);

		}
	}
}

void euler_field::Boundary_Macroscopic()
{
	//left
	for (int y = 1; y < Ny ; y++)
	{
		for (int k = 0; k < Q; k++)
		{
			double rho = euler_node[1][y]->Get_rho();
			euler_node[0][y]->Set_rho(rho);

			//boundary
			//euler_node[0][y]->Set_Velocity(_u_inlet,0);
			//boundary

			double feq = FeqStandardD2Q9(k, euler_node[0][y]->Get_rho(), euler_node[0][y]->Get_Velocity_x(), euler_node[0][y]->Get_Velocity_y());

			double fneq1 = euler_node[1][y]->Get_New_f(k);
			double fneq2 = FeqStandardD2Q9(k, euler_node[1][y]->Get_rho(), euler_node[1][y]->Get_Velocity_x(), euler_node[1][y]->Get_Velocity_y());

			double result = feq + (fneq1 - fneq2)*(1 - 1/_tau);


			euler_node[0][y]->Set_f(result,k);

		}
	}
	//right
	for (int y = 1; y < Ny; y++)
	{
		for (int k = 0; k < Q; k++)
		{

			double rho = euler_node[Nx-1][y]->Get_rho();
			euler_node[Nx][y]->Set_rho(rho);

			//boundary
			//euler_node[Nx][y]->Set_Velocity(_u_inlet, 0);
			//boundary



			double feq = FeqStandardD2Q9(k, euler_node[Nx][y]->Get_rho(), euler_node[Nx][y]->Get_Velocity_x(), euler_node[Nx][y]->Get_Velocity_y());


			double fneq1 = euler_node[Nx-1][y]->Get_New_f(k);
			double fneq2 = FeqStandardD2Q9(k, euler_node[Nx - 1][y]->Get_rho(), euler_node[Nx - 1][y]->Get_Velocity_x(), euler_node[Nx - 1][y]->Get_Velocity_y());

			double result = feq + (fneq1 - fneq2) *(1 - 1 / _tau);


			euler_node[Nx][y]->Set_f(result, k);


		}
	}

	//top
	for (int x = 0; x <= Nx; x++)
	{
		for (int k = 0; k < Q; k++)
		{
			double rho = euler_node[x][Ny - 1]->Get_rho();
			euler_node[x][Ny]->Set_rho(rho);



			//boundary
			//euler_node[x][Ny]->Set_Velocity(_u_inlet, 0);
			//boundary
			double feq  = FeqStandardD2Q9(k, euler_node[x][Ny]->Get_rho(), euler_node[x][Ny]->Get_Velocity_x(), euler_node[x][Ny]->Get_Velocity_y());
			
			double fneq1 = euler_node[x][Ny - 1]->Get_f(k);
			double fneq2 = FeqStandardD2Q9(k, euler_node[x][Ny - 1]->Get_rho(), euler_node[x][Ny - 1]->Get_Velocity_x(), euler_node[x][Ny - 1]->Get_Velocity_y());

			double result = feq + (fneq1 - fneq2) *(1 - 1 / _tau);

			euler_node[x][Ny]->Set_f(result, k);
		}
	}



	//bottom

	for (int x = 0; x <= Nx; x++)
	{
		for (int k = 0; k < Q; k++)
		{
			double rho = euler_node[x][1]->Get_rho();
			euler_node[x][0]->Set_rho(rho);
			//boundary
			//euler_node[x][0]->Set_Velocity(_u_inlet, 0);
			//boundary
			double feq = FeqStandardD2Q9(k, euler_node[x][0]->Get_rho(), euler_node[x][0]->Get_Velocity_x(), euler_node[x][0]->Get_Velocity_y());

			double fneq1 = euler_node[x][1]->Get_f(k);
			double fneq2 = FeqStandardD2Q9(k, euler_node[x][1]->Get_rho(), euler_node[x][1]->Get_Velocity_x(), euler_node[x][1]->Get_Velocity_y());

			double result = feq + (fneq1 - fneq2)  * (1 - 1 / _tau);

			euler_node[x][0]->Set_f(result, k);

		}
		//cout << euler_node[x][1]->Get_Velocity_x()<<endl;
	}
	//cout << endl << endl;

}

void euler_field::Boundary_Condition()
{

	//left
	for (int y = 1; y < Ny; y++)
	{
		for (int k = 0; k < Q; k++)
		{

			//boundary
			euler_node[0][y]->Set_Velocity(_u_inlet,0);
			//boundary


		}
	}
	//right
	for (int y = 1; y < Ny; y++)
	{
		for (int k = 0; k < Q; k++)
		{
			euler_node[Nx][y]->Set_Velocity(_u_inlet, 0);
		}
	}

	//top
	for (int x = 0; x <= Nx; x++)
	{
		for (int k = 0; k < Q; k++)
		{
			euler_node[x][Ny]->Set_Velocity(0, 0);
		}
	}



	//bottom

	for (int x = 0; x <= Nx; x++)
	{
		for (int k = 0; k < Q; k++)
		{
			euler_node[x][0]->Set_Velocity(0 , 0);
		}
		//cout << euler_node[x][1]->Get_Velocity_x()<<endl;
	}




}

void euler_field::Top_Boundary_Velocity_Set(double u, double v)
{
	for (int x = 0; x <= Nx; x++)
	{
		for (int k = 0; k < Q; k++)
		{
			euler_node[x][Ny]->Set_Velocity(u, v);
		}
	}
}

void euler_field::Button_Velocity_Boundary_Set(double u, double v)
{
	for (int x = 0; x <= Nx; x++)
	{
		for (int k = 0; k < Q; k++)
		{
			euler_node[x][0]->Set_Velocity(u, 0);
		}
		//cout << euler_node[x][1]->Get_Velocity_x()<<endl;
	}

}

void euler_field::left_Velocity_Boundary_Set(double u, double v)
{
	for (int y = 1; y < Ny; y++)
	{
		for (int k = 0; k < Q; k++)
		{

			//boundary
			euler_node[0][y]->Set_Velocity(u, v);
			//boundary


		}
	}
}

void euler_field::right_Velocity_Boundary_Set(double u, double v)
{
	for (int y = 1; y < Ny; y++)
	{
		for (int k = 0; k < Q; k++)
		{
			euler_node[Nx][y]->Set_Velocity(u, v);
		}
	}
}


void euler_field::Output_grid()
{
	ostringstream name2;
	name2 << "grid" << ".plt";
	ofstream out(name2.str().c_str());

	out << "filetype = grid " << endl;
	out << "variables =" << " \"x \" " << "," << " \" y \" " << endl;
	out << "zone i = " << Nx +1<< "\t" << "j = " << Ny+1;
	out << "\t" << "f = point";

	out << endl;

	for (int i = 0; i <= Nx; i++)
	{
		for (int j = 0; j <= Ny; j++)
		{

			out << euler_node[i][j]->Get_Position_x();
			out << "\t";
			out << euler_node[i][j]->Get_Position_y();


			out << endl;
		}

	}

}





void euler_field::Output_field(double value, int step)
{
	ostringstream name2;
	name2 << "Re=" << value << "_euler_node_" << step << ".plt";
	ofstream out(name2.str().c_str());



	out << "filetype = solution" << endl;
	out << "variables =" << " \"u\" " << "," << " \"v \" " << "," << " \"P \" " << endl;
	out << "zone i = " << Nx+1 << "\t" << "j = " << Ny+1 << "\t" << "f = point";

	out << endl;

	for (int i = 0; i <= Nx; i++)
	{
		for (int j = 0; j <= Ny; j++)
		{


			out << euler_node[i][j]->Get_Velocity_x() << "\t";
			out << euler_node[i][j]->Get_Velocity_y() << "\t";
			out << euler_node[i][j]->Get_rho() << "\t";
			//out << euler_node[i][j]->Get_Dentisy() << "\t";

			//out << euler_node[i][j]->Get_u_noF() << "\t";
			//out << euler_node[i][j]->Get_v_noF() << "\t";

			//out << euler_node[i][j]->Get_Body_force_fx() << "\t";
			//out << euler_node[i][j]->Get_Body_force_fy() << "\t";

			for (int m = 0; m < _NF; m++)
			{

			}


			out << endl;
		}

	}

}

void euler_field::Output_Fluid_NF(int NF_m, int NF_step, int cylinder_step)
{
	ostringstream name;


	//name << "m = " << NF_m << "_lagrange_node" << cylinder_step <<"forcing_step="<< NF_step <<".dat";

	name << "euler_node" << cylinder_step << "m = " << NF_m << "forcing_step=" << NF_step << ".csv";


	ofstream out(name.str().c_str());



	//out << euler_node[31][54]->Get_NF_u(1) << "\t";
	//out << euler_node[31][54]->Get_NF_v(1) << "\t";
	//out << euler_node[31][54]->Get_NF_Fx(1) << "\t";
	//out << euler_node[31][54]->Get_NF_Fy(1) << "\t";
	//out << endl;

	for (int i = 0; i <= Nx; i++)
	{
		for (int j = 0; j <= Ny; j++)
		{


			out << euler_node[i][j]->Get_Position_x() << "\t";
			out << euler_node[i][j]->Get_Position_y() << "\t";
			out << endl;

			for (int m = 0; m < _NF; m++)
			{
				out << "m = " << m << "\t";
				out <<euler_node[i][j]->Get_NF_u(m) << "\t";
				out << euler_node[i][j]->Get_NF_v(m) << "\t";
				out << euler_node[i][j]->Get_NF_Fx(m) << "\t";
				out << euler_node[i][j]->Get_NF_Fy(m) << "\t";
				out << endl;
			}



			out << "result" << endl;


			out << euler_node[i][j]->Get_NF_u(_NF-1) << "\t";
			out << euler_node[i][j]->Get_NF_v(_NF -1) << "\t";

			out << euler_node[i][j]->Get_Velocity_x() << "\t";
			out << euler_node[i][j]->Get_Velocity_y() << "\t";

			out << euler_node[i][j]->Get_Body_force_fx() << "\t";
			out << euler_node[i][j]->Get_Body_force_fy() << "\t";
			out << endl;

		}


		
	}


}

void euler_field::Output_all(double value, int step)
{
	ostringstream name2;
	name2 << "Re=" << value << "_euler_node_" << step << "all.plt";
	ofstream out(name2.str().c_str());



	out << "filetype = solution" << endl;
	out << "variables =" << " \"u\" " << "," << " \"v \" " << "," << " \"P \" " << endl;
	out << "zone i = " << Nx + 1 << "\t" << "j = " << Ny + 1 << "\t" << "f = point";

	out << endl;

	for (int i = 0; i <= Nx; i++)
	{
		for (int j = 0; j <= Ny; j++)
		{
			out << euler_node[i][j]->Get_Position_x() << "\t";
			out << euler_node[i][j]->Get_Position_y() << "\t";

			out << euler_node[i][j]->Get_Velocity_x() << "\t";
			out << euler_node[i][j]->Get_Velocity_y() << "\t";

			out << euler_node[i][j]->Get_Body_force_fx() << "\t";
			out << euler_node[i][j]->Get_Body_force_fy() << "\t";


			//out << euler_node[i][j]->Get_rho() << "\t";
			//out << euler_node[i][j]->Get_Dentisy() << "\t";

			//out << euler_node[i][j]->Get_u_noF() << "\t";
			//out << euler_node[i][j]->Get_v_noF() << "\t";

			//out << euler_node[i][j]->Get_Body_force_fx() << "\t";
			//out << euler_node[i][j]->Get_Body_force_fy() << "\t";

			for (int m = 0; m < _NF; m++)
			{

			}


			out << endl;
		}

	}


}

void euler_field::Clear_NF_F_ij()
{
	for (int i = 0; i < euler_node.size(); i++)
	{
		for (int j = 0; j < euler_node[0].size(); j++)
		{
			for (int m = 0; i < _NF; i++)
			{
				euler_node[i][j]->Set_NF_F(m, 0, 0);
			}
		}
	}


}

void euler_field::Clear_NF_u()
{
	for (int i = 0; i < euler_node.size(); i++)
	{
		for (int j = 0; j < euler_node[0].size(); j++)
		{
			for (int m = 0; m < _NF; m++)
			{
				euler_node[i][j]->Set_NF_u(m, 0, 0);
				

			}

		}
	}



}

void euler_field::Output_f(int x, int y)
{
	cout << "("<< x <<","<< y << ")"<< endl;
	cout << euler_node[x][y]->Get_f(6) << "\t";
	cout << euler_node[x][y]->Get_f(2) << "\t";
	cout << euler_node[x][y]->Get_f(5) << "\t";
	cout << endl;
	cout << euler_node[x][y]->Get_f(3) << "\t";
	cout << euler_node[x][y]->Get_f(0) << "\t";
	cout << euler_node[x][y]->Get_f(1) << "\t";
	cout << endl;
	cout << euler_node[x][y]->Get_f(7) << "\t";
	cout << euler_node[x][y]->Get_f(4) << "\t";
	cout << euler_node[x][y]->Get_f(8) << "\t";
	cout << endl;
	cout << "__________________________________________________________________________________" << endl;



}

void euler_field::Output_Parameter()
{
	cout << "Nx=" << Nx << endl;
	cout << "Ny=" << Ny << endl;
	cout << "L=" << _L << endl;
	cout << "Re=" << _Re << endl;
	cout << "tau=" << _tau << endl;
	cout << "u_inlet=" << _u_inlet << endl;
}





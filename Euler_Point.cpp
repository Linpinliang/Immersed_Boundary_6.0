#include "Euler_Point.h"


Euler_Point::Euler_Point()
{
}

Euler_Point::Euler_Point(int x, int y)
{

	Set_Position(x, y);

	Set_Velocity(0, 0);

	Set_rho(1.0);

	Set_Temperature(1.0);

	Set_Body_force(0, 0);



	//初始化是计算分布函数
	Collision();
	Set_u_v_noF(0, 0);
}


Euler_Point::~Euler_Point()
{
}

void Euler_Point::macroscopic()
{
	//cout << "密度:" << Rho << endl;
	//计算密度
	Rho0 = Rho;
	Rho = 0;
	for (int k = 0; k < Q; k++) {
		Rho += f[k];
		// << "密度:" << Rho << endl;
	}


	u0 = u;
	v0 = v;


	//计算速度
	Set_Velocity(0, 0);
	for (int k = 0; k < Q; k++) {
		u += e[k][0] * f[k];
		v += e[k][1] * f[k];
		//	Print_Velocity();
	}

	u /= Rho;
	v /= Rho;

	//Print_Velocity();



}

void Euler_Point::Collision()
{
	for (int k = 0; k < Q; k++) {
		f[k] = FeqStandardD2Q9(k, Rho, u, v);
		//cout << f[k] << endl;
	}
}








void Euler_Point::Set_Velocity(double vx, double vy)
{
	this->u = vx;
	this->v = vy;
}

void Euler_Point::Set_Velocity0(double u, double v)
{
	u0 = u;
	v0 = v;
}

void Euler_Point::Set_Body_force(double fx, double fy)
{
	this->forcing_x = fx;
	this->forcing_y = fy;
}

void Euler_Point::Set_Position(int x, int y)
{
	this->grid_x = x;
	this->grid_y = y;
}

void Euler_Point::Set_Temperature(double t)
{
	this->T = t;
}

void Euler_Point::Set_Pressure(double p)
{

}

void Euler_Point::Set_rho(double rho)
{
	this->Rho = rho;
}

void Euler_Point::Set_u_v_noF(double u, double v)
{
	u_noF = u;
	v_noF = v;
}

void Euler_Point::Set_f(double value, int N)
{
	f[N] = value;
}

void Euler_Point::Set_New_f(double value, int N)
{
	New_f[N] = value;
}

void Euler_Point::Set_F(double value, int N)
{
	F[N] = value;
}






double Euler_Point::Get_Body_force_fx()
{
	return forcing_x;
}

double Euler_Point::Get_Body_force_fy()
{
	return forcing_y;
}

double Euler_Point::Get_Velocity_x()
{
	return u;
}

double Euler_Point::Get_Velocity_y()
{
	return v;
}

double Euler_Point::Get_Velocity0_x()
{
	return u0;
}

double Euler_Point::Get_Velocity0_y()
{
	return v0;
}

int Euler_Point::Get_Position_x()
{
	return grid_x;
}

int Euler_Point::Get_Position_y()
{
	return grid_y;
}

double Euler_Point::Get_u_noF()
{
	return u_noF;
}

double Euler_Point::Get_v_noF()
{
	return v_noF;
}

double Euler_Point::Get_Pressure()
{
	return Rho / 3;
}

double Euler_Point::Get_Temperature()
{
	return T;
}

double Euler_Point::Get_rho()
{
	return Rho;
}

double Euler_Point::Get_f(int N)
{
	return f[N];
}

double Euler_Point::Get_New_f(int N)
{
	return New_f[N];
}

double Euler_Point::Get_F(int N)
{
	return F[N];
}



void Euler_Point::Print_f()
{
	cout << "Print_f()" << endl;
	cout << f[6] << "\t" << f[2] << "\t" << f[5] << endl;
	cout << f[3] << "\t" << f[0] << "\t" << f[1] << endl;
	cout << f[7] << "\t" << f[4] << "\t" << f[8] << endl;
}

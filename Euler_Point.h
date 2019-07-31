#ifndef EULER_POINT_H
#define EULER_POINT_H
#include "Parameter.h"
#include <iostream>
#include "Function.h"
/*


*/

using namespace std;

class Euler_Point
{
private:
	//网格
	int grid_x;
	int grid_y;



	//速度
	double u;
	double v;

	double u0;
	double v0;

	//unfoced velocity 
	double u_noF;
	double v_noF;


	//网格类型
	char style;

	//温度
	double T;
	double T0;
	//密度
	double Rho;
	double Rho0;

	double P;

	//体力项
	double forcing_x;
	double forcing_y;

	//分布函数
	double New_f[9];
	double f[9];

	//力项
	double F[9];


public:

	Euler_Point();

	//初始化
	Euler_Point(int x, int y);

	~Euler_Point();





	//LBM

	//计算微观量
	void macroscopic();
	//碰撞
	void Collision();
	//力项
	//void Scond_Forcing_Step();

	//void Collision_With_Body_Force();





	//Set函数
	//设置速度
	void Set_Velocity(double vx, double vy);

	void Set_Velocity0(double u ,double v);
	//设置体积力项;
	void Set_Body_force(double fx, double fy);
	//设置点的位置
	void Set_Position(int x, int y);
	//设置温度
	void Set_Temperature(double t);
	//设置压力
	void Set_Pressure(double p);
	//设置密度
	void Set_rho(double rho);
	void Set_u_v_noF(double u, double v);




	void Set_f(double value, int N);
	void Set_New_f(double value, int N);
	void Set_F(double value, int N);






	//Get函数

	double Get_Body_force_fx();
	double Get_Body_force_fy();

	double Get_Velocity_x();
	double Get_Velocity_y();

	double Get_Velocity0_x();
	double Get_Velocity0_y(); 


	int Get_Position_x();
	int Get_Position_y();

	double Get_u_noF();
	double Get_v_noF();


	double Get_Pressure();
	double Get_Temperature();
	double Get_rho();

	double Get_f(int N);

	double  Get_New_f(int N);
	double Get_F(int N);



	//Print函数
	//输出分布函数
	void Print_f();




};


#endif // !EULER_POINT_H

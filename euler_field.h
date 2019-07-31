#ifndef EULER_FIELD_H
#define  EULER_FIELD_H
#include <vector>
#include "Euler_Point.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <omp.h>


class euler_field
{
public:
	int Nx;
	int Ny;
	int _Max_Step;
	vector<vector<Euler_Point*> > euler_node;


	double _u_inlet;
	//雷诺数
	double _Re;
	//松弛因子
	double _tau;
	//特征尺度
	double _L;


public:
	void Set_Re(double Re);
	void Set_L(double L);
	void Set_tau();
	void Set_u_inlet(double u);
	void Set_Max_step(int N);


	double Get_Re();
	double Get_L();
	double Get_tau();
	double Get_u_inlet();
	int Get_Max_set();

	//euler_field();

	euler_field();
	euler_field(int x,int y,double Re ,double L ,double u_inlet,int _max_step);
	~euler_field();


	void evolution();

	void Collision();

	void Stream();

	void Macroscopic();

	void Boundary_Macroscopic();

	void Boundary_Condition();

	void Top_Boundary_Velocity_Set(double u ,double v);

	void Button_Velocity_Boundary_Set(double u, double v);

	void left_Velocity_Boundary_Set(double u, double v);

	void right_Velocity_Boundary_Set(double u, double v);

	
	//Output 
	void Output_grid();

	void Output_field(double value, int step);
	
	void Output_f(int x, int y);

	void Output_Parameter();

	void output_all(double value, int step);
};


#endif // EULER_FIELD_H


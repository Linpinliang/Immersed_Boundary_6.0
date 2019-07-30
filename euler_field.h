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

	int _output_step;


	double _u_inlet;
	//雷诺数
	double _Re;
	//松弛因子
	double _tau;
	//特征尺度
	double _L;
	//mluti
	int _NF;


public:
	void Set_Re(double Re);
	void Set_L(double L);
	void Set_tau();
	void Set_u_inlet(double u);
	void Set_Max_step(int N);
	void Set_NF(int nf);


	double Get_Re();
	double Get_L();
	double Get_tau();
	double Get_u_inlet();
	int Get_Max_step();
	int Get_NF();



	//euler_field();

	euler_field();
	euler_field(int x, int y, double Re, double L, double u_inlet, int _max_step);

	euler_field(int x, int y, double Re, double L, double u_inlet, int _max_step,int nf_number);


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

	void Output_Fluid_NF(int NF_m, int NF_step, int cylinder_step);

	void Output_all(double value, int step);


	void Clear_NF_F_ij();
	void Clear_NF_u();


	
	void Output_f(int x, int y);

	void Output_Parameter();
};


#endif // EULER_FIELD_H


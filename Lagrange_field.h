#ifndef LAGRANGE_FIELD_H
#define LAGRANGE_FIELD_H
#include <vector>
#include "Lagrange_Point.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "Function.h"
#define pi  3.14159265358979323846



using namespace std;

class Lagrange_field
{
public:
	double _center_x;
	double _center_y;
	double _r;
	int _Number_of_lagrange_node;
	double _delta_Sb;
	vector<Lagrange_Point*> lagrange_node;

	//multi method parameter

	int _NF;


public:
	//
	Lagrange_field();
	
	//Ô²Öù
	Lagrange_field(double center_x,double center_y, double r,int number_of_node);


	//multi method constructor
	Lagrange_field(double center_x, double center_y, double r, int number_of_node,int NF);


	void Set_center_x(double x);
	void Set_center_y(double y);
	void Set_Number_of_lagrange_node(int N);
	void Set_delta_Sb();
	void Set_r(double r);
	void Set_lagrange_node();
	
	//multi 
	void Set_NF_lagrange_node(int nf_number);

	//clear NF parameter
	void Clear_NF_Fb();
	void Clear_NF_ub();


	double Get_center_x(double x);
	double Get_center_y(double y);
	int Get_Number_of_lagrange_node(int N);
	double Get_r(double r);
	double Get_delta_Sb();

	//Output Data
	void Output_Solid(double value, int step);
	void Output_Solid_NF(int NF_m, int NF_step ,int cylinder_step);

	~Lagrange_field();



};

#endif 
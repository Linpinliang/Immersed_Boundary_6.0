#ifndef LAGRANGE_POINT_H
#define LAGRANGE_POINT_H
#include <vector>
#include <iostream>
using namespace std;


class Lagrange_Point
{
private:
	double _x;
	double _y;

	double _ub_noF;
	double _vb_noF;

	double _Ub;
	double _Vb;

	double _Fx;
	double _Fy;



	//NF method parameter
	int _NF;

	vector<double > _NF_Fx;
	vector<double > _Nf_Fy;
	vector<double > _NF_ub;
	vector<double > _NF_vb;


public:
	//Lagrange_Point();
	Lagrange_Point(double x, double y);

	//multi_mothod_constructor
	Lagrange_Point(double x, double y ,int NF);






	//Set method
	void Set_Position(double x, double y);
	void Set_ub_vb_noF(double ub_noF, double vb_noF);
	void Set_Ub_Vb(double Ub, double Vb);
	void Set_F(double Fx, double Fy);


	//multi_method Parameter
	void Set_NF(int nf);
	void Set_NF_ub_vb(double ub ,double vb,int m);
	void Set_NF_F(double Fx ,double Fy ,int m);

	//Get method
	double Get_Position_x();
	double Get_Position_y();

	double get_ub_nof();
	double get_vb_nof();

	double Get_Ub();
	double Get_Vb();

	double Get_Fx();
	double Get_Fy();

	//multi_method get
	double Get_NF_ub(int m);
	double Get_NF_vb(int m);

	double Get_NF_Fx(int m);
	double Get_NF_Fy(int m);

	int Get_NF_size();







	~Lagrange_Point();
};


#endif //  

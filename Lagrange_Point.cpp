#include "Lagrange_Point.h"


/*
Lagrange_Point::Lagrange_Point()
{
}
*/
Lagrange_Point::Lagrange_Point(double x, double y)
{
	Set_Position(x, y);
	Set_ub_vb_noF(0, 0);
	Set_Ub_Vb(0, 0);
	Set_F(0, 0);
}

void Lagrange_Point::Set_Position(double x, double y)
{
	_x = x;
	_y = y;
}

void Lagrange_Point::Set_ub_vb_noF(double ub_noF, double vb_noF)
{
	_ub_noF = ub_noF;
	_vb_noF = vb_noF;
}

void Lagrange_Point::Set_Ub_Vb(double Ub, double Vb)
{
	_Ub = Ub;
	_Vb = Vb;

}

void Lagrange_Point::Set_F(double Fx, double Fy)
{
	_Fx = Fx;
	_Fy = Fy;
}

void Lagrange_Point::Set_deleta_sb(double DELTA_SB)
{
	_delta_sb = DELTA_SB;
}



double Lagrange_Point::Get_Position_x()
{
	return _x;
}

double Lagrange_Point::Get_Position_y()
{
	return _y;
}

double Lagrange_Point::get_ub_nof()
{
	return _ub_noF;
}

double Lagrange_Point::get_vb_nof()
{
	return _vb_noF;
}

double Lagrange_Point::Get_Ub()
{
	return _Ub;
}

double Lagrange_Point::Get_Vb()
{
	return _Vb;
}

double Lagrange_Point::Get_Fx()
{
	return _Fx;
}

double Lagrange_Point::Get_Fy()
{
	return _Fy;
}

double Lagrange_Point::Get_delta_sb()
{
	return _delta_sb;
}


Lagrange_Point::~Lagrange_Point()
{
}

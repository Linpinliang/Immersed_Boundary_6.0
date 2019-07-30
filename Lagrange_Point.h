#ifndef LAGRANGE_POINT_H
#define LAGRANGE_POINT_H
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

	double _delta_sb;

	

public:
	//Lagrange_Point();
	Lagrange_Point(double x, double y);




	void Set_Position(double x, double y);
	void Set_ub_vb_noF(double ub_noF, double vb_noF);
	void Set_Ub_Vb(double Ub, double Vb);
	void Set_F(double Fx, double Fy);
	void Set_deleta_sb(double DELTA_SB);


	double Get_Position_x();
	double Get_Position_y();

	double get_ub_nof();
	double get_vb_nof();

	double Get_Ub();
	double Get_Vb();

	double Get_Fx();
	double Get_Fy();
	double Get_delta_sb();
	

	~Lagrange_Point();
};


#endif //  

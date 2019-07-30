#include "Function.h"



Function::Function()
{
}


Function::~Function()
{
}





double FeqStandardD2Q9(int k, double rho, double Vx, double Vy)
{
	double eu, uv, feq;
	eu = (e[k][0] * Vx + e[k][1] * Vy);
	uv = (Vx *  Vx + Vy * Vy);
	feq = w[k] * rho*(1.0 + 3.0*eu + 4.5*eu*eu - 1.5*uv);

	return feq;
}

double d_h(double r) {
	double dh;
	if (fabs(r) < 1) {

		dh = 3 - 2 * fabs(r) + sqrt(1 + 4 * fabs(r) - 4 * r * r);
		dh = dh / 8;
		//cout << "case1"<<endl;
		return dh;

	}
	if (fabs(r) < 2 && fabs(r) >= 1) {

		dh = 5 - 2 * fabs(r) - sqrt(-7 + 12 * fabs(r) - 4 * r*r);

		dh = dh / 8;
		//cout << 2<<endl;
		return dh;
	}
	if (fabs(r) >= 0) {
		//cout << 3<<endl;
		return 0;
	}



	return 0;

}


double force_distribution_function(double Fx, double Fy, double u, double v, int k)
{
	/*
	double result1 = 3 * w[k] * (e[k][0] * Fx + e[k][1] * Fy - (u*Fx) - (v*Fy));
	double result2 = 9 * w[k] * ((e[k][0] * u + e[k][1] * v) + (e[k][0] * Fx + e[k][1] * Fy));

	return result1 + result2;
	*/


	double result = 3 * w[k] * (e[k][0] * Fx + e[k][1] * Fy);
	return result;
}


double D_function(double x1, double y1, double x2, double y2) {


	double h = 1;

	double result = 0;
	result = d_h((x1 - x2) / h) * d_h((y1 - y2) / h) / h / h;

	return result;
}



double Distance(double x1, double y1, double x2, double y2) {

	return sqrt((x1 - x2) * (x1 - x2) + (y1 - y2)*(y1 - y2));
}


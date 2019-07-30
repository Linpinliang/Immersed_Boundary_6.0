#ifndef FEQSTANDARDD2Q9_H
#define FEQSTANDARDD2Q9_H
#include "Parameter.h"
#include <cmath>
class Function
{
public:
	Function();
	~Function();
};


double FeqStandardD2Q9(int k, double rho, double Vx, double Vy);

double force_distribution_function(double Fx, double Fy, double u, double v, int k);

double d_h(double r);

double D_function(double x1, double y1, double x2, double y2);

double Distance(double x1, double y1, double x2, double y2);

#endif




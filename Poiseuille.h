#ifndef POISEUILLE_H
#define POISEUILLE_H



#include "euler_field.h"

void Poiseuille01() {
	int grid_x = 10;
	int grid_y = 20;
	int Re = 1;
	int L = double(grid_y);

	euler_field flow1 = euler_field(grid_x, grid_y, Re, L, 0.1, 1001);

	flow1.Top_Boundary_Velocity_Set(0, 0);
	flow1.Button_Velocity_Boundary_Set(0, 0);
	flow1.left_Velocity_Boundary_Set(0.1, 0);
	flow1.right_Velocity_Boundary_Set(0.1, 0);

	flow1.evolution();

}




#endif // !POISEUILLE_H

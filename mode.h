
#ifndef MODE_H
#define MODE_H
#include"tokamak.h"

typedef struct MODE
{
	int n;
	int m;
	int pa,pb;
	double r_s;
	double delta_r;
}Mode;

void G_R_theta(Grid * const grid, Tokamak * const tok,Mode * const pmode,double **G_2D);

#endif

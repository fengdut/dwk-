#ifndef DIST_H
#define DIST_H

void F0_3D(const Slowing* slow,const Grid *grid, double *** F_3D, double ***FE_3D, double ***FR_3D);
void Lambda_b_L_3D(const Grid *grid, const Tokamak *tok,double*** lambda_b_3D,double *** b_lambda_3D);
void Theta(double ***b_lambda_3D,const Grid *grid, double ***Theta_3D);
#endif


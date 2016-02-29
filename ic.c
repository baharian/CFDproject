/*
 *  ATMS 502 / CSE 566 -- Fall, 2012
 *  pgm5:  3D nonlinear quasi-compressible flow
 *  SOHEIL BAHARIAN
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "defs.h"

/* ============================ ic =========================================
 * IC sets the initial condition
 *
 * Arguments:
 *	t1, t2, p1, p2		real arrays		past (n = -1) and current (n = 0)
 *										temperature and pressure fields
 *	u1, v1, w1			real arrays		past (n = -1) flow fields
 *	u2, v2, w2			real arrays		current (n = 0) flow fields
 *	rho_t, rho_w		real arrays		densities at $theta$ and $w$ levels
 *	dx, dy, dz			real			grid spacing
 *	i1, i2x, i2y, i2z	integers		indices bounding array data
 *	i2x_u, i2y_v, i2z_w	integers		upper bounds of velocity field data
 *	nx, ny, nz			integers		number of grid points
 *	uPrt				real			max random U perturbation along X
 *	tPrt				real array		temperature perturbations info
 *	tBar, g, Cp			reals			physical constants
 */
void ic(t1, t2, p1, p2, u1, v1, w1, u2, v2, w2, rho_t, rho_w, dx, dy, dz, i1, i2x, i2y, i2z, i2x_u, i2y_v, i2z_w, nx, ny, nz, uPrt, tPrt, tBar, g, Cp)
float	t1[NXDIM][NYDIM][NZDIM], t2[NXDIM][NYDIM][NZDIM], p1[NXDIM][NYDIM][NZDIM], p2[NXDIM][NYDIM][NZDIM],
		u1[NXDIM_U][NYDIM][NZDIM], v1[NXDIM][NYDIM_V][NZDIM], w1[NXDIM][NYDIM][NZDIM_W],
		u2[NXDIM_U][NYDIM][NZDIM], v2[NXDIM][NYDIM_V][NZDIM], w2[NXDIM][NYDIM][NZDIM_W],
		rho_t[NZDIM], rho_w[NZDIM_W], uPrt, tPrt[prt][8], dx, dy, dz, tBar, g, Cp;
int		i1, i2x, i2y, i2z, i2x_u, i2y_v, i2z_w, nx, ny, nz;
{
	float pi = 4.0 * atan(1.0);
	float R = 287.0, P0 = powf(10.0, 5.0);

	int i, j, k, m;
	float x, y, z, r, x0, y0, z0, rx, ry, rz, T, P;
	
	/* base state */
	for (k = i1; k <= i2z; k++)
	{
		z = dz * ((float) (k - i1) + 0.5);	/* height at $theta$-level */
		T = tBar - (g/Cp) * z;
		P = P0 * powf(T/tBar, Cp/R);
		rho_t[k] = P/(R*T);	/* density at $theta$-level */

		if (k > i1)
			rho_w[k] = (rho_t[k] + rho_t[k-1]) / 2.0;	/* density at $w$-level */
	}
	rho_w[i1] = 0.0;	/* density at z = 0, corresponding to k = 1 */
	rho_w[i2z_w] = 0.0;	/* density at z = NZ*dz, corresponding to k = NZ+1 */

	#pragma omp parallel for private(i, j, k, m, x, y, z, x0, y0, z0, rx, ry, rz, r) shared(t1, p1, t2, p2, tPrt)
	for (i = i1; i <= i2x; i++)
		for (j = i1; j <= i2y; j++)
			for (k = i1; k <= i2z; k++)
			{
				t1[i][j][k] = tBar;
				p1[i][j][k] = 0.0;

				x = dx * ((float) (i - i1) + 0.5);
				y = dy * ((float) (j - i1) + 0.5);
				z = dz * ((float) (k - i1) + 0.5);

				for (m = 0; m < prt; m++)
				{
					/* perturbation center */
					x0 = tPrt[m][2];
					y0 = tPrt[m][3];
					z0 = tPrt[m][4];
					/* perturbation radii */
					rx = tPrt[m][5];
					ry = tPrt[m][6];
					rz = tPrt[m][7];

					r = sqrt(powf((x - x0) / rx, 2.0) + powf((y - y0) / ry, 2.0) + powf((z - z0) / rz, 2.0));

					if (r <= 1.0)
						t1[i][j][k] += tPrt[m][0] * (cos(pi * r) + 1.0) / 2.0;
				}

				t2[i][j][k] = t1[i][j][k];
				p2[i][j][k] = p1[i][j][k];
			}

	/* x-component of the wind field, u(x,y,z), is random */
	srand(0.0);

	#pragma omp parallel for private(i, j, k) shared(u1, u2)
	for (i = i1; i <= i2x_u; i++)
		for (j = i1; j <= i2y; j++)
			for (k = i1; k <= i2z; k++)
			{
				u1[i][j][k] = uPrt * ((float) rand() / (float) RAND_MAX - 0.5);
				u2[i][j][k] = u1[i][j][k];
			}

	/* y-component of the wind field, v(x,y,z) */
	#pragma omp parallel for private(i, j, k, m, x, y, z, x0, y0, z0, rx, ry, rz, r) shared(v1, v2, tPrt)
	for (i = i1; i <= i2x; i++)
		for (j = i1; j <= i2y_v; j++)
			for (k = i1; k <= i2z; k++)
			{
				v1[i][j][k] = 0.0;

				x = dx * ((float) (i - i1) + 0.5);
				y = dy * ((float) (j - i1) + 0.5);
				z = dz * ((float) (k - i1) + 0.5);

				for (m = 0; m < prt; m++)
				{
					/* perturbation center */
					x0 = tPrt[m][2];
					y0 = tPrt[m][3];
					z0 = tPrt[m][4];
					/* perturbation radii */
					rx = tPrt[m][5];
					ry = tPrt[m][6];
					rz = tPrt[m][7];

					r = sqrt(powf((x - x0) / rx, 2.0) + powf((y - y0) / ry, 2.0) + powf((z - z0) / rz, 2.0));

					if (r <= 1.0)
						v1[i][j][k] += tPrt[m][1] * (cos(pi * r) + 1.0) / 2.0;
				}

				v2[i][j][k] = v1[i][j][k];
			}

	/* z-component of the wind field, w(x,y,z), is zero initially */
	#pragma omp parallel for private(i, j, k) shared(w1, w2)
	for (i = i1; i <= i2x; i++)
		for (j = i1; j <= i2y; j++)
			for (k = i1; k <= i2z_w; k++)
			{
				w1[i][j][k] = 0.0;
				w2[i][j][k] = w1[i][j][k];
			}

	return;
}
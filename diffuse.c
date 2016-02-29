/*
 *  ATMS 502 / CSE 566 -- Fall, 2012
 *  pgm5:  3D nonlinear quasi-compressible flow
 *  SOHEIL BAHARIAN
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "defs.h"

/* ======================= diffuse =========================================
 * DIFFUSE calculates the diffusive effects on the fields
 *
 * Arguments:
 *	t2, t3				real arrays		values at current and next steps
 *	u1, v1, w1			real arrays		past X, Y, Z velocity flow fields
 *	u3, v3, w3			real arrays		future X, Y, Z velocity flow fields
 *	K, Kt				reals			diffusion coefficients for
 *										velocity and temperature fields
 *	dx, dy, dz			real			grid spacing in x and y directions
 *	dt, tstep			real			time steps
 *	i1, i2x, i2y, i2z	integers		indices bounding array data
 *	i2x_u, i2y_v, i2z_w	integers		upper bounds of velocity field data
 *	nx, ny, nz			integers		number of grid points
 */
void diffuse(t2, t3, u1, v1, w1, u3, v3, w3, K, Kt, dx, dy, dz, dt, tstep, i1, i2x, i2y, i2z, i2x_u, i2y_v, i2z_w, nx, ny, nz)
float	t2[NXDIM][NYDIM][NZDIM], t3[NXDIM][NYDIM][NZDIM],
		u1[NXDIM_U][NYDIM][NZDIM], v1[NXDIM][NYDIM_V][NZDIM], w1[NXDIM][NYDIM][NZDIM_W],
		u3[NXDIM_U][NYDIM][NZDIM], v3[NXDIM][NYDIM_V][NZDIM], w3[NXDIM][NYDIM][NZDIM_W],
		K, Kt, dx, dy, dz, dt, tstep;
int		i1, i2x, i2y, i2z, i2x_u, i2y_v, i2z_w, nx, ny, nz;
{
	int i, j, k;

	/* diffusion effects for u(x,y,z) */
	#pragma omp parallel for private(i, j, k) shared(u1, u3)
	for (i = i1 + 1; i <= i2x_u - 1; i++)	/* NOTE: boundary conditions set u[I1] = - u[I1+1] and u[I2X_U] = - u[I2X_U-1] */
		for (j = i1; j <= i2y; j++)
			for (k = i1; k <= i2z; k++)
				u3[i][j][k] += tstep * K
					* ((u1[i+1][j][k] - 2.0 * u1[i][j][k] + u1[i-1][j][k]) / (dx * dx)
					 + (u1[i][j+1][k] - 2.0 * u1[i][j][k] + u1[i][j-1][k]) / (dy * dy)
					 + (u1[i][j][k+1] - 2.0 * u1[i][j][k] + u1[i][j][k-1]) / (dz * dz));

	/* diffusion effects for v(x,y,z) */
	#pragma omp parallel for private(i, j, k) shared(v1, v3)
	for (i = i1; i <= i2x; i++)
		for (j = i1; j <= i2y_v - 1; j++)	/* NOTE: boundary conditions set v[I2Y_V] = v[I1] */
			for (k = i1; k <= i2z; k++)
				v3[i][j][k] += tstep * K
					* ((v1[i+1][j][k] - 2.0 * v1[i][j][k] + v1[i-1][j][k]) / (dx * dx)
					 + (v1[i][j+1][k] - 2.0 * v1[i][j][k] + v1[i][j-1][k]) / (dy * dy)
					 + (v1[i][j][k+1] - 2.0 * v1[i][j][k] + v1[i][j][k-1]) / (dz * dz));

	/* diffusion effects for w(x,y,z) */
	#pragma omp parallel for private(i, j, k) shared(w1, w3)
	for (i = i1; i <= i2x; i++)
		for (j = i1; j <= i2y; j++)
			for (k = i1 + 1; k <= i2z_w - 1; k++)	/* NOTE: boundary conditions set w[i1] = w[i2z_w] = 0 at all times */
				w3[i][j][k] += tstep * K
					* ((w1[i+1][j][k] - 2.0 * w1[i][j][k] + w1[i-1][j][k]) / (dx * dx)
					 + (w1[i][j+1][k] - 2.0 * w1[i][j][k] + w1[i][j-1][k]) / (dy * dy)
					 + (w1[i][j][k+1] - 2.0 * w1[i][j][k] + w1[i][j][k-1]) / (dz * dz));

	/* diffusion effects for temperature */
	#pragma omp parallel for private(i, j, k) shared(t2, t3)
	for (i = i1; i <= i2x; i++)
		for (j = i1; j <= i2y; j++)
			for (k = i1; k <= i2z; k++)
				t3[i][j][k] += dt * Kt
					* ((t2[i+1][j][k] - 2.0 * t2[i][j][k] + t2[i-1][j][k]) / (dx * dx)
					 + (t2[i][j+1][k] - 2.0 * t2[i][j][k] + t2[i][j-1][k]) / (dy * dy)
					 + (t2[i][j][k+1] - 2.0 * t2[i][j][k] + t2[i][j][k-1]) / (dz * dz));

	return;
}
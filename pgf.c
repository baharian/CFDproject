/*
 *  ATMS 502 / CSE 566 -- Fall, 2012
 *  pgm5:  3D nonlinear quasi-compressible flow
 *  SOHEIL BAHARIAN
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "defs.h"

/* ======================= PGF =============================================
 * PGF calculates the pressure gradiant forces and the bouyancy effect
 *
 * Arguments:
 *	t2					real array		temperature at current step
 *	p1, p3				real arrays		pressure at previous and next steps
 *	u3, v3, w3			real arrays		future X, Y, Z velocity flow fields
 *	rho_t, rho_w		real arrays		densities at $theta$ and $w$ levels
 *	dx, dy, dz			real			grid spacing in x and y directions
 *	tstep					real			time step
 *	i1, i2x, i2y, i2z	integers		indices bounding array data
 *	i2x_u, i2y_v, i2z_w	integers		upper bounds of velocity field data
 *	nx, ny, nz			integers		number of grid points
 *	tBar, g, Cs			reals			physical constants
 */
void pgf(t2, p1, p3, u3, v3, w3, rho_t, rho_w, dx, dy, dz, tstep, i1, i2x, i2y, i2z, i2x_u, i2y_v, i2z_w, nx, ny, nz, tBar, g, Cs)
float	t2[NXDIM][NYDIM][NZDIM], p1[NXDIM][NYDIM][NZDIM], p3[NXDIM][NYDIM][NZDIM],
		u3[NXDIM_U][NYDIM][NZDIM], v3[NXDIM][NYDIM_V][NZDIM], w3[NXDIM][NYDIM][NZDIM_W],
		rho_t[NZDIM], rho_w[NZDIM_W], dx, dy, dz, tstep, tBar, g, Cs;
int		i1, i2x, i2y, i2z, i2x_u, i2y_v, i2z_w, nx, ny, nz;
{
	int i, j, k;

	/* pressure-gradiant effects for u(x,y,z) */
	#pragma omp parallel for private(i, j, k) shared(p1, u3, rho_t)
	for (i = i1 + 1; i <= i2x_u - 1; i++)	/* NOTE: boundary conditions set u[I1] = - u[I1+1] and u[I2X_U] = - u[I2X_U-1] */
		for (j = i1; j <= i2y; j++)
			for (k = i1; k <= i2z; k++)
				u3[i][j][k] += - tstep * ((p1[i][j][k] - p1[i-1][j][k]) / dx) / rho_t[k];

	/* pressure-gradiant effects for v(x,y,z) */
	#pragma omp parallel for private(i, j, k) shared(p1, v3, rho_t)
	for (i = i1; i <= i2x; i++)
		for (j = i1; j <= i2y_v - 1; j++)	/* NOTE: boundary conditions set v[I2Y_V] = v[I1] */
			for (k = i1; k <= i2z; k++)
				v3[i][j][k] += - tstep * ((p1[i][j][k] - p1[i][j-1][k]) / dy) / rho_t[k];

	/* pressure-gradiant and bouyancy effects for w(x,y,z) */
	#pragma omp parallel for private(i, j, k) shared(p1, t2, w3, rho_w)
	for (i = i1; i <= i2x; i++)
		for (j = i1; j <= i2y; j++)
			for (k = i1 + 1; k <= i2z_w - 1; k++)	/* NOTE: boundary conditions set w[i1] = w[i2z_w] = 0 at all times */
				w3[i][j][k] += tstep
					* (	- ((p1[i][j][k] - p1[i][j][k-1]) / dz) / rho_w[k]	/* rho_w[i1] = rho_w[i2z_w] = 0, but we never touch these! */
						+ (g / tBar) * ((t2[i][j][k] - tBar) + (t2[i][j][k-1] - tBar)) / 2.0);

	/* this is to set the correct boundary conditions for u3, v3, and w3. note that t2
	   (current step) and p1 (previous step) have the correct boundaries already. */
	bc(t2, p1, u3, v3, w3, i1, i2x, i2y, i2z, i2x_u, i2y_v, i2z_w, nx, ny, nz);

	/* time-evolution for the pressure field, p(x,y,z) */
	#pragma omp parallel for private(i, j, k) shared(p1, u3, v3, w3, p3, rho_t, rho_w)
	for (i = i1; i <= i2x; i++)
		for (j = i1; j <= i2y; j++)
			for (k = i1; k <= i2z; k++)
				p3[i][j][k] = p1[i][j][k] - tstep * (Cs * Cs)
					* (rho_t[k] * (u3[i+1][j][k] - u3[i][j][k]) / dx + rho_t[k] * (v3[i][j+1][k] - v3[i][j][k]) / dy + (rho_w[k+1] * w3[i][j][k+1] - rho_w[k] * w3[i][j][k]) / dz);

	return;
}
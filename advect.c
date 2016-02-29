/*
 *  ATMS 502 / CSE 566 -- Fall, 2012
 *  pgm5:  3D nonlinear quasi-compressible flow
 *  SOHEIL BAHARIAN
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "defs.h"

/* ======================= advect ===========================================
 * ADVECT moves the fields forward by one time step
 *
 * Arguments:
 *	t2, t3				real arrays		values at current and next steps
 *	u1, v1, w1			real arrays		past X, Y, Z velocity flow fields
 *	u2, v2, w2			real arrays		current X, Y, Z velocity flow fields
 *	u3, v3, w3			real arrays		future X, Y, Z velocity flow fields
 *	dx, dy, dz			real			grid spacing in x and y directions
 *	dt, tstep			real			time step
 *	i1, i2x, i2y, i2z	integers		indices bounding array data
 *	i2x_u, i2y_v, i2z_w	integers		upper bounds of velocity field data
 *	nx, ny, nz			integers		number of grid points
 */
void advect(t2, t3, u1, v1, w1, u2, v2, w2, u3, v3, w3, dx, dy, dz, dt, tstep, i1, i2x, i2y, i2z, i2x_u, i2y_v, i2z_w, nx, ny, nz)
float	t2[NXDIM][NYDIM][NZDIM], t3[NXDIM][NYDIM][NZDIM],
		u1[NXDIM_U][NYDIM][NZDIM], v1[NXDIM][NYDIM_V][NZDIM], w1[NXDIM][NYDIM][NZDIM_W],
		u2[NXDIM_U][NYDIM][NZDIM], v2[NXDIM][NYDIM_V][NZDIM], w2[NXDIM][NYDIM][NZDIM_W],
		u3[NXDIM_U][NYDIM][NZDIM], v3[NXDIM][NYDIM_V][NZDIM], w3[NXDIM][NYDIM][NZDIM_W],
		dx, dy, dz, dt, tstep;
int		i1, i2x, i2y, i2z, i2x_u, i2y_v, i2z_w, nx, ny, nz;
{
	void advect1D();

	int i, j, k;
	float temp1[NXDIM], temp2[NXDIM], vel1D[NX+1];	/* CAREFUL! check the size of vel1D */

	/* advect the x-component of the wind field, u(x,y,z) */
	#pragma omp parallel for private(i, j, k) shared(u1, u2, v2, w2, u3)
	for (i = i1 + 1; i <= i2x_u - 1; i++)	/* NOTE: boundary conditions set u[I1] = - u[I1+1] and u[I2X_U] = - u[I2X_U-1] */
		for (j = i1; j <= i2y; j++)
			for (k = i1; k <= i2z; k++)
				u3[i][j][k] = u1[i][j][k] - (tstep / 4.0)
					* ((u2[i+1][j][k] * u2[i+1][j][k] - u2[i-1][j][k] * u2[i-1][j][k]) / dx
					 + ((v2[i][j+1][k] + v2[i-1][j+1][k]) * (u2[i][j+1][k] - u2[i][j][k]) + (v2[i][j][k] + v2[i-1][j][k]) * (u2[i][j][k] - u2[i][j-1][k])) / dy
					 + ((w2[i][j][k+1] + w2[i-1][j][k+1]) * (u2[i][j][k+1] - u2[i][j][k]) + (w2[i][j][k] + w2[i-1][j][k]) * (u2[i][j][k] - u2[i][j][k-1])) / dz);
	
	/* advect the y-component of the wind field, v(x,y,z) */
	#pragma omp parallel for private(i, j, k) shared(v1, u2, v2, w2, v3)
	for (i = i1; i <= i2x; i++)
		for (j = i1; j <= i2y_v - 1; j++) /* NOTE: boundary conditions set v[I2Y_V] = v[I1] */
			for (k = i1; k <= i2z; k++)
				v3[i][j][k] = v1[i][j][k] - (tstep / 4.0)
					* (((u2[i+1][j][k] + u2[i+1][j-1][k]) * (v2[i+1][j][k] - v2[i][j][k]) + (u2[i][j][k] + u2[i][j-1][k]) * (v2[i][j][k] - v2[i-1][j][k])) / dx
					 + (v2[i][j+1][k] * v2[i][j+1][k] - v2[i][j-1][k] * v2[i][j-1][k]) / dy
					 + ((w2[i][j][k+1] + w2[i][j-1][k+1]) * (v2[i][j][k+1] - v2[i][j][k]) + (w2[i][j][k] + w2[i][j-1][k]) * (v2[i][j][k] - v2[i][j][k-1])) / dz);

	/* advect the z-component of the wind field, w(x,y,z) */
	#pragma omp parallel for private(i, j, k) shared(w1, u2, v2, w2, w3)
	for (i = i1; i <= i2x; i++)
		for (j = i1; j <= i2y; j++)
			for (k = i1 + 1; k <= i2z_w - 1; k++)	/* NOTE: boundary conditions set w[I1] = w[I2Z_W] = 0 at all times */
				w3[i][j][k] = w1[i][j][k] - (tstep / 4.0)
					* (((u2[i+1][j][k] + u2[i+1][j][k-1]) * (w2[i+1][j][k] - w2[i][j][k]) + (u2[i][j][k] + u2[i][j][k-1]) * (w2[i][j][k] - w2[i-1][j][k])) / dx
					 + ((v2[i][j+1][k] + v2[i][j+1][k-1]) * (w2[i][j+1][k] - w2[i][j][k]) + (v2[i][j][k] + v2[i][j][k-1]) * (w2[i][j][k] - w2[i][j-1][k])) / dy
					 + (w2[i][j][k+1] * w2[i][j][k+1] - w2[i][j][k-1] * w2[i][j][k-1]) / dz);

	/* advect the temperature field, t(x,y,z)... */
	/* ...advect by (dt/2) along X direction... */
	for (j = i1; j <= i2y; j++)
		for (k = i1; k <= i2z; k++)
		{
			#pragma omp parallel for private(i) shared(temp1, temp2, t2)
			for (i = 0; i < NXDIM; i++)
			{
				temp1[i] = t2[i][j][k];
				temp2[i] = 0.0;	/* initialize the second array */
			}
			#pragma omp parallel for private(i) shared(vel1D, u2)
			for (i = i1; i <= i2x_u; i++)
				vel1D[i - i1] = u2[i][j][k];

			advect1D(temp1, temp2, vel1D, dx, dt / 2.0, i1, i2x, nx, 'P');	/* 'P' is Piecewise-Linear */

			#pragma omp parallel for private(i) shared(temp2, t3)
			for (i = i1; i <= i2x; i++)
				t3[i][j][k] = temp2[i];
		}

	/* set the correct boundary conditions for t3 */
	/* NOTE: there is no pressure field in this routine, so t3 is used again for the second
			 argument of bc(); and u1, v1, w1 are also dummy variables passed to bc() */
	bc(t3, t3, u1, v1, w1, i1, i2x, i2y, i2z, i2x_u, i2y_v, i2z_w, nx, ny, nz);

	/* ...advect by (dt/2) along Y direction... */
	for (i = i1; i <= i2x; i++)
		for (k = i1; k <= i2z; k++)
		{
			#pragma omp parallel for private(j) shared(temp1, temp2, t3)
			for (j = 0; j < NYDIM; j++)
			{
				temp1[j] = t3[i][j][k];
				temp2[j] = 0.0;	/* initialize the second array */
			}
			#pragma omp parallel for private(j) shared(vel1D, v2)
			for (j = i1; j <= i2y_v; j++)
				vel1D[j - i1] = v2[i][j][k];

			advect1D(temp1, temp2, vel1D, dy, dt / 2.0, i1, i2y, ny, 'P');	/* 'P' is Piecewise-Linear */

			#pragma omp parallel for private(j) shared(temp2, t3)
			for (j = i1; j <= i2y; j++)
				t3[i][j][k] = temp2[j];
		}

	/* set the correct boundary conditions for t3 */
	/* NOTE: there is no pressure field in this routine, so t3 is used again for the second
			 argument of bc(); and u1, v1, w1 are also dummy variables passed to bc() */
	bc(t3, t3, u1, v1, w1, i1, i2x, i2y, i2z, i2x_u, i2y_v, i2z_w, nx, ny, nz);

	/* ...advect by (dt) along Z direction... */
	for (i = i1; i <= i2x; i++)
		for (j = i1; j <= i2y; j++)
		{
			#pragma omp parallel for private(k) shared(temp1, temp2, t3)
			for (k = 0; k < NZDIM; k++)
			{
				temp1[k] = t3[i][j][k];
				temp2[k] = 0.0;	/* initialize the second array */
			}
			#pragma omp parallel for private(k) shared(vel1D, w2)
			for (k = i1; k <= i2z_w; k++)
				vel1D[k - i1] = w2[i][j][k];

			advect1D(temp1, temp2, vel1D, dz, dt, i1, i2z, nz, 'P');	/* 'P' is Piecewise-Linear */

			#pragma omp parallel for private(k) shared(temp2, t3)
			for (k = i1; k <= i2z; k++)
				t3[i][j][k] = temp2[k];
		}

	/* set the correct boundary conditions for t3 */
	/* NOTE: there is no pressure field in this routine, so t3 is used again for the second
			 argument of bc(); and u1, v1, w1 are also dummy variables passed to bc() */
	bc(t3, t3, u1, v1, w1, i1, i2x, i2y, i2z, i2x_u, i2y_v, i2z_w, nx, ny, nz);

	/* ...advect by (dt/2) along Y direction... */
	for (i = i1; i <= i2x; i++)
		for (k = i1; k <= i2z; k++)
		{
			#pragma omp parallel for private(j) shared(temp1, temp2, t3)
			for (j = 0; j < NYDIM; j++)
			{
				temp1[j] = t3[i][j][k];
				temp2[j] = 0.0;	/* initialize the second array */
			}
			#pragma omp parallel for private(j) shared(vel1D, v2)
			for (j = i1; j <= i2y_v; j++)
				vel1D[j - i1] = v2[i][j][k];

			advect1D(temp1, temp2, vel1D, dy, dt / 2.0, i1, i2y, ny, 'P');	/* 'P' is Piecewise-Linear */

			#pragma omp parallel for private(j) shared(temp2, t3)
			for (j = i1; j <= i2y; j++)
				t3[i][j][k] = temp2[j];
		}

	/* set the correct boundary conditions for t3 */
	/* NOTE: there is no pressure field in this routine, so t3 is used again for the second
			 argument of bc(); and u1, v1, w1 are also dummy variables passed to bc() */
	bc(t3, t3, u1, v1, w1, i1, i2x, i2y, i2z, i2x_u, i2y_v, i2z_w, nx, ny, nz);

	/* ...advect by (dt/2) along X direction */
	for (j = i1; j <= i2y; j++)
		for (k = i1; k <= i2z; k++)
		{
			#pragma omp parallel for private(i) shared(temp1, temp2, t3)
			for (i = 0; i < NXDIM; i++)
			{
				temp1[i] = t3[i][j][k];
				temp2[i] = 0.0;	/* initialize the second array */
			}
			#pragma omp parallel for private(i) shared(vel1D, u2)
			for (i = i1; i <= i2x_u; i++)
				vel1D[i - i1] = u2[i][j][k];

			advect1D(temp1, temp2, vel1D, dx, dt / 2.0, i1, i2x, nx, 'P');	/* 'P' is Piecewise-Linear */

			#pragma omp parallel for private(i) shared(temp2, t3)
			for (i = i1; i <= i2x; i++)
				t3[i][j][k] = temp2[i];
		}

	return;
}
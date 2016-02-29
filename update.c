/*
 *  ATMS 502 / CSE 566 -- Fall, 2012
 *  pgm5:  3D nonlinear quasi-compressible flow
 *  SOHEIL BAHARIAN
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "defs.h"

/* ========================= update ===================================================
 * UPDATE replaces old values with new ones
 *
 * Arguments:
 *	t1, t2, t3			real array		temperature at past, current, and future steps
 *	p1, p2, p3			real arrays		pressure at past, current, and future steps
 *	u1, v1, w1			real arrays		past X, Y, Z velocity flow fields
 *	u2, v2, w2			real arrays		current X, Y, Z velocity flow fields
 *	u3, v3, w3			real arrays		future X, Y, Z velocity flow fields
 *	i1, i2x, i2y, i2z	integers		indices bounding array data
 *	i2x_u, i2y_v, i2z_w	integers		upper bounds of velocity field data
 *	nx, ny, nz			integers		number of grid points
 *	n					integer			current time step
 */
void update(t1, t2, t3, p1, p2, p3, u1, v1, w1, u2, v2, w2, u3, v3, w3, i1, i2x, i2y, i2z, i2x_u, i2y_v, i2z_w, nx, ny, nz, n)
float	t1[NXDIM][NYDIM][NZDIM], t2[NXDIM][NYDIM][NZDIM], t3[NXDIM][NYDIM][NZDIM],
		p1[NXDIM][NYDIM][NZDIM], p2[NXDIM][NYDIM][NZDIM], p3[NXDIM][NYDIM][NZDIM],
		u1[NXDIM_U][NYDIM][NZDIM], v1[NXDIM][NYDIM_V][NZDIM], w1[NXDIM][NYDIM][NZDIM_W],
		u2[NXDIM_U][NYDIM][NZDIM], v2[NXDIM][NYDIM_V][NZDIM], w2[NXDIM][NYDIM][NZDIM_W],
		u3[NXDIM_U][NYDIM][NZDIM], v3[NXDIM][NYDIM_V][NZDIM], w3[NXDIM][NYDIM][NZDIM_W];
int		i1, i2x, i2y, i2z, i2x_u, i2y_v, i2z_w, nx, ny, nz, n;
{
	int i, j, k;

	if (n > 1)
	{
		#pragma omp parallel for private(i, j, k) shared(t1, p1, t2, p2)
		for (i = i1; i <= i2x; i++)
			for (j = i1; j <= i2y; j++)
				for (k = i1; k <= i2z; k++)
				{
					t1[i][j][k] = t2[i][j][k];
					p1[i][j][k] = p2[i][j][k];
				}

		#pragma omp parallel for private(i, j, k) shared(u1, u2)
		for (i = i1; i <= i2x_u; i++)
			for (j = i1; j <= i2y; j++)
				for (k = i1; k <= i2z; k++)
					u1[i][j][k] = u2[i][j][k];

		#pragma omp parallel for private(i, j, k) shared(v1, v2)
		for (i = i1; i <= i2x; i++)
			for (j = i1; j <= i2y_v; j++)
				for (k = i1; k <= i2z; k++)
					v1[i][j][k] = v2[i][j][k];

		#pragma omp parallel for private(i, j, k) shared(w1, w2)
		for (i = i1; i <= i2x; i++)
			for (j = i1; j <= i2y; j++)
				for (k = i1; k <= i2z_w; k++)
					w1[i][j][k] = w2[i][j][k];
	}

	#pragma omp parallel for private(i, j, k) shared(t2, p2, t3, p3)
	for (i = i1; i <= i2x; i++)
		for (j = i1; j <= i2y; j++)
			for (k = i1; k <= i2z; k++)
			{
				t2[i][j][k] = t3[i][j][k];
				p2[i][j][k] = p3[i][j][k];
			}

	#pragma omp parallel for private(i, j, k) shared(u2, u3)
	for (i = i1; i <= i2x_u; i++)
		for (j = i1; j <= i2y; j++)
			for (k = i1; k <= i2z; k++)
				u2[i][j][k] = u3[i][j][k];

	#pragma omp parallel for private(i, j, k) shared(v2, v3)
	for (i = i1; i <= i2x; i++)
		for (j = i1; j <= i2y_v; j++)
			for (k = i1; k <= i2z; k++)
				v2[i][j][k] = v3[i][j][k];

	#pragma omp parallel for private(i, j, k) shared(w2, w3)
	for (i = i1; i <= i2x; i++)
		for (j = i1; j <= i2y; j++)
			for (k = i1; k <= i2z_w; k++)
				w2[i][j][k] = w3[i][j][k];

	return;
}
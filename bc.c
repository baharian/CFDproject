/*
 *  ATMS 502 / CSE 566 -- Fall, 2012
 *  pgm5:  3D nonlinear quasi-compressible flow
 *  SOHEIL BAHARIAN
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "defs.h"

/* ============================ bc =========================================
 * BC sets the boundary conditions
 *
 * Arguments:
 *	t2, p2				real arrays		current (n = 1) temperature
 *										and pressure fields
 *	u2, v2, w2			real arrays		current (n = 1) flow fields
 *	i1, i2x, i2y, i2z	integers		indices bounding array data
 *	i2x_u, i2y_v, i2z_w	integers		upper bounds of velocity field data
 *	nx, ny, nz			integers		number of grid points
 */
void bc(t2, p2, u2, v2, w2, i1, i2x, i2y, i2z, i2x_u, i2y_v, i2z_w, nx, ny, nz)
float	t2[NXDIM][NYDIM][NZDIM], p2[NXDIM][NYDIM][NZDIM],
		u2[NXDIM_U][NYDIM][NZDIM], v2[NXDIM][NYDIM_V][NZDIM], w2[NXDIM][NYDIM][NZDIM_W];
int		i1, i2x, i2y, i2z, i2x_u, i2y_v, i2z_w, nx, ny, nz;
{
	/* NOTE: this routine is written for any (general) number of ghost points */
	int i, j, k;

	/* rigid lids for w at k = I1 and at k = I2Z_W */
	#pragma omp parallel for private(i, j) shared(w2)
	for (i = i1; i <= i2x; i++)
		for (j = i1; j <= i2y; j++)
		{
			w2[i][j][i1] = 0.0;
			w2[i][j][i2z_w] = 0.0;
		}

	/* left and right (X) */
	for (i = 1; i <= BC_WIDTH; i++)
	{
		#pragma omp parallel for private(j, k) shared(t2, p2, u2, v2, w2)
		for (j = i1; j <= i2y; j++)
			for (k = i1; k <= i2z; k++)
			{
				/* u is antisymmetric in x-direction */
				/* CAREFUL! first and last elements, u[0] and u[NXDIM_U-1], are not set by this method... */
				u2[i1 - i + 1][j][k] = - u2[i1 + i][j][k];
				u2[i2x_u + i - 1][j][k] = - u2[i2x_u - i][j][k];
				/* ...they happen for i = BC_WIDTH+1, and now they're set! (but i don't think it's necessary) */
				u2[0][j][k] = - u2[2 * i1 + 1][j][k];
				u2[NXDIM_U-1][j][k] = - u2[i2x_u - i1 - 1][j][k];
				
				/* v, w, t, p are symmetric in x-direction */
				v2[i1 - i][j][k] = v2[i1 + i][j][k];
				v2[i2x + i][j][k] = v2[i2x - i][j][k];

				w2[i1 - i][j][k] = w2[i1 + i][j][k];
				w2[i2x + i][j][k] = w2[i2x - i][j][k];

				t2[i1 - i][j][k] = t2[i1 + i][j][k];
				t2[i2x + i][j][k] = t2[i2x - i][j][k];
				
				p2[i1 - i][j][k] = p2[i1 + i][j][k];
				p2[i2x + i][j][k] = p2[i2x - i][j][k];
			}

		/* set the back-most v's */
		#pragma omp parallel for private(k) shared(v2)
		for (k = i1; k <= i2z; k++)
		{
			v2[i1 - i][i2y_v][k] = v2[i1 + i][i2y_v][k];
			v2[i2x + i][i2y_v][k] = v2[i2x - i][i2y_v][k];
		}

		/* set the top-most w's */
		#pragma omp parallel for private(j) shared(w2)
		for (j = i1; j <= i2y; j++)
		{
			w2[i1 - i][j][i2z_w] = w2[i1 + i][j][i2z_w];
			w2[i2x + i][j][i2z_w] = w2[i2x - i][j][i2z_w];
		}
	}

	/* front and back (Y) */
	for (j = 1; j <= BC_WIDTH; j++)
	{
		#pragma omp parallel for private(i, k) shared(t2, p2, u2, v2, w2)
		for (i = i1; i <= i2x; i++)
			for (k = i1; k <= i2z; k++)
			{
				/* u, v, w, t, p are periodic in y-direction */
				u2[i][i1 - j][k] = u2[i][i2y + 1 - j][k];
				u2[i][i2y + j][k] = u2[i][i1 - 1 + j][k];

				v2[i][i1 - j][k] = v2[i][i2y_v - j][k];
				v2[i][i2y_v + j][k] = v2[i][i1 + j][k];
				v2[i][i2y_v][k] = v2[i][i1][k];	/* the boundary is at j = I1 and j = I2Y_V */

				w2[i][i1 - j][k] = w2[i][i2y + 1 - j][k];
				w2[i][i2y + j][k] = w2[i][i1 - 1 + j][k];

				t2[i][i1 - j][k] = t2[i][i2y + 1 - j][k];
				t2[i][i2y + j][k] = t2[i][i1 - 1 + j][k];

				p2[i][i1 - j][k] = p2[i][i2y + 1 - j][k];
				p2[i][i2y + j][k] = p2[i][i1 - 1 + j][k];
			}

		/* set the right-most u's */
		#pragma omp parallel for private(k) shared(u2)
		for (k = i1; k <= i2z; k++)
		{
			u2[i2x_u][i1 - j][k] = u2[i2x_u][i2y + 1 - j][k];
			u2[i2x_u][i2y + j][k] = u2[i2x_u][i1 - 1 + j][k];
		}

		/* set the top-most w's */
		#pragma omp parallel for private(i) shared(w2)
		for (i = i1; i <= i2x; i++)
		{
			w2[i][i1 - j][i2z_w] = w2[i][i2y + 1 - j][i2z_w];
			w2[i][i2y + j][i2z_w] = w2[i][i1 - 1 + j][i2z_w];
		}
	}

	/* bottom and top (Z) */
	for (k = 1; k <= BC_WIDTH; k++)
	{
		#pragma omp parallel for private(i, j) shared(t2, p2, u2, v2, w2)
		for (i = i1; i <= i2x; i++)
			for (j = i1; j <= i2y; j++)
			{
				/* u, v, w, t, p are zero-gradient in z-direction */
				u2[i][j][i1 - k] = u2[i][j][i1];
				u2[i][j][i2z + k] = u2[i][j][i2z];

				v2[i][j][i1 - k] = v2[i][j][i1];
				v2[i][j][i2z + k] = v2[i][j][i2z];

				w2[i][j][i1 - k] = w2[i][j][i1];
				w2[i][j][i2z_w + k] = w2[i][j][i2z_w];

				t2[i][j][i1 - k] = t2[i][j][i1];
				t2[i][j][i2z + k] = t2[i][j][i2z];

				p2[i][j][i1 - k] = p2[i][j][i1];
				p2[i][j][i2z + k] = p2[i][j][i2z];
			}

		/* set the right-most u's */
		#pragma omp parallel for private(j) shared(u2)
		for (j = i1; j <= i2y; j++)
		{
			u2[i2x_u][j][i1 - k] = u2[i2x_u][j][i1];
			u2[i2x_u][j][i2z + k] = u2[i2x_u][j][i2z];
		}

		/* set the back-most v's */
		#pragma omp parallel for private(i) shared(v2)
		for (i = i1; i <= i2x; i++)
		{
			v2[i][i2y_v][i1 - k] = v2[i][i2y_v][i1];
			v2[i][i2y_v][i2z + k] = v2[i][i2y_v][i2z];
		}
	}

	return;
}
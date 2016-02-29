/*
 *  ATMS 502 / CSE 566 -- Fall, 2012
 *  pgm5:  3D nonlinear quasi-compressible flow
 *  SOHEIL BAHARIAN
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "defs.h"

/* ========================= stats ==========================================
 * STATS computes and prints out the max and min values
 *
 * Arguments:
 *	t2, p2				real arrays		current (n = 1) temperature
 *										and pressure fields
 *	u2, v2, w2			real arrays		current (n = 1) flow fields
 *	i1, i2x, i2y, i2z	integers		indices bounding array data
 *	i2x_u, i2y_v, i2z_w	integers		upper bounds of velocity field data
 *	nx, ny, nz			integers		number of grid points
 *	n					integer			time step counter
 *	dt					real			time step value
 *	tBar				real			ambient temperature
 */
void stats(t2, p2, u2, v2, w2, i1, i2x, i2y, i2z, i2x_u, i2y_v, i2z_w, nx, ny, nz, n, dt, tBar)
float	t2[NXDIM][NYDIM][NZDIM], p2[NXDIM][NYDIM][NZDIM],
		u2[NXDIM_U][NYDIM][NZDIM], v2[NXDIM][NYDIM_V][NZDIM], w2[NXDIM][NYDIM][NZDIM_W];
float	dt, tBar;
int		i1, i2x, i2y, i2z, i2x_u, i2y_v, i2z_w, nx, ny, nz, n;
{
	int i, j, k;
	float	tMinTemp = t2[i1][i1][i1] - tBar, tMaxTemp = t2[i1][i1][i1] - tBar, pMinTemp = p2[i1][i1][i1], pMaxTemp = p2[i1][i1][i1],
			uMinTemp = u2[i1][i1][i1], uMaxTemp = u2[i1][i1][i1], vMinTemp = v2[i1][i1][i1], vMaxTemp = v2[i1][i1][i1], wMinTemp = w2[i1][i1][i1], wMaxTemp = w2[i1][i1][i1];
	
	for (i = i1; i <= i2x; i++)
		for (j = i1; j <= i2y; j++)
			for (k = i1; k <= i2z; k++)
			{
				/* find tmax */
				if ((t2[i][j][k] - tBar) > tMaxTemp)
					tMaxTemp = t2[i][j][k] - tBar;

				/* find tmin */
				if ((t2[i][j][k] - tBar) < tMinTemp)
					tMinTemp = t2[i][j][k] - tBar;

				/* find pmax */
				if (p2[i][j][k] > pMaxTemp)
					pMaxTemp = p2[i][j][k];

				/* find pmin */
				if (p2[i][j][k] < pMinTemp)
					pMinTemp = p2[i][j][k];
			}

	for (i = i1 + 1; i <= i2x_u - 1; i++)	/* NOTE: boundary conditions set u[I1] = - u[I1+1] and u[I2X_U] = - u[I2X_U-1] */
		for (j = i1; j <= i2y; j++)
			for (k = i1; k <= i2z; k++)
			{
				/* find umax */
				if (u2[i][j][k] > uMaxTemp)
					uMaxTemp = u2[i][j][k];

				/* find umin */
				if (u2[i][j][k] < uMinTemp)
					uMinTemp = u2[i][j][k];
			}

	for (i = i1; i <= i2x; i++)
		for (j = i1; j <= i2y_v - 1; j++)	/* NOTE: boundary conditions set v[I2Y_V] = v[I1] */
			for (k = i1; k <= i2z; k++)
			{
				/* find vmax */
				if (v2[i][j][k] > vMaxTemp)
					vMaxTemp = v2[i][j][k];

				/* find vmin */
				if (v2[i][j][k] < vMinTemp)
					vMinTemp = v2[i][j][k];
			}

	for (i = i1; i <= i2x; i++)
		for (j = i1; j <= i2y; j++)
			for (k = i1 + 1; k <= i2z_w - 1; k++)	/* NOTE: boundary conditions set w[i1] = w[i2z_w] = 0 at all times */
			{
				/* find wmax */
				if (w2[i][j][k] > wMaxTemp)
					wMaxTemp = w2[i][j][k];

				/* find wmin */
				if (w2[i][j][k] < wMinTemp)
					wMinTemp = w2[i][j][k];
			}

	printf("%4d  %3.2f  %3.3f  %3.3f  %3.3f  %3.3f  %3.3f  %3.3f  %3.4f  %3.4f  %4.2f  %4.2f\n"
		, n, dt * (float) n, uMinTemp, uMaxTemp, vMinTemp, vMaxTemp, wMinTemp, wMaxTemp, tMinTemp, tMaxTemp, pMinTemp, pMaxTemp);

	return;
}
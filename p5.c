/*
 *  ATMS 502 / CSE 566 -- Fall, 2012
 *  pgm5:  3D nonlinear quasi-compressible flow
 *  SOHEIL BAHARIAN
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "defs.h"

int main(void)
{
	char *name = "SOHEIL BAHARIAN";

	/* CAREFUL! check with online cases */
	float dx = 50.0;
	float dy = 50.0;
	float dz = 50.0;
	float dt = 0.15;

	/* arrays and other variables */
	float	t1[NXDIM][NYDIM][NZDIM], t2[NXDIM][NYDIM][NZDIM], t3[NXDIM][NYDIM][NZDIM],
			p1[NXDIM][NYDIM][NZDIM], p2[NXDIM][NYDIM][NZDIM], p3[NXDIM][NYDIM][NZDIM],
			u1[NXDIM_U][NYDIM][NZDIM], v1[NXDIM][NYDIM_V][NZDIM], w1[NXDIM][NYDIM][NZDIM_W],
			u2[NXDIM_U][NYDIM][NZDIM], v2[NXDIM][NYDIM_V][NZDIM], w2[NXDIM][NYDIM][NZDIM_W],
			u3[NXDIM_U][NYDIM][NZDIM], v3[NXDIM][NYDIM_V][NZDIM], w3[NXDIM][NYDIM][NZDIM_W],
			rho_t[NZDIM], rho_w[NZDIM_W],
			tPlot[NX][NY][NZ], pPlot[NX][NY][NZ], uPlot[NX][NY][NZ], vPlot[NX][NY][NZ], wPlot[NX][NY][NZ];
	int		i, j, k, n, nstep, nplot;
	float	K, Kt, Cs, tBar = 300.0, g = 9.81, Cp = 1004.0, tstep;
	float	uPrt, tPrt[prt][8];	/* center @ (X,Y,Z) + radius in (X,Y,Z) directions + T-value + V-value = 8 variables */

	/* definitions for routines */
	void putfield();
	void ic(), bc(), advect(), diffuse(), pgf(), update(), stats();

	/* read-in parameters */
	/*for (n = 1; n <= prt; n++)
	{
		printf("Perturbation %d, temperature value: ", n);
		scanf("%f", &(tPrt[n-1][0]));

		printf("Perturbation %d, Y-velocity value: ", n);
		scanf("%f", &(tPrt[n-1][1]));

		printf("Perturbation %d, center (X,Y,Z) in meters: ", n);
		scanf("%f %f %f", &(tPrt[n-1][2]), &(tPrt[n-1][3]), &(tPrt[n-1][4]));

		printf("Perturbation %d, radius (rX,rY,rZ) in meters: ", n);
		scanf("%f %f %f", &(tPrt[n-1][5]), &(tPrt[n-1][6]), &(tPrt[n-1][7]));
	}

	printf("Enter the diffusion coefficient for the velocity fields: ");
	scanf("%f", &K);
	printf("Enter the diffusion coefficient for the temperature field: ");
	scanf("%f", &Kt);

	printf("Enter the speed of sound: ");
	scanf("%f", &Cs);

	printf("Enter number of time steps to take: ");
	scanf("%d", &nstep);
	printf("Enter plot interval, in steps: ");
	scanf("%d", &nplot);*/

	/* CAREFUL! check with online cases */
	tPrt[0][0] = -15.0;	/* T perturbation */
	tPrt[0][1] = -35.0;	/* V perturbation */
	tPrt[0][2] = 25.0;		tPrt[0][3] = 7775.0;	tPrt[0][4] = 2050.0;	/* center @ (X,Y,Z) */
	tPrt[0][5] = 3500.0;	tPrt[0][6] = 999999.0;	tPrt[0][7] = 3500.0;	/* radius in (X,Y,Z) */

	tPrt[1][0] = -25.0;	/* T perturbation */
	tPrt[1][1] = 35.0;	/* V perturbation */
	tPrt[1][2] = 11225.0;	tPrt[1][3] = 4225.0;	tPrt[1][4] = 4400.0;	/* center @ (X,Y,Z) */
	tPrt[1][5] = 2000.0;	tPrt[1][6] = 2000.0;	tPrt[1][7] = 2000.0;	/* radius in (X,Y,Z) */

	uPrt = 1.0;

	K = 15.0;
	Kt = 2.0;

	Cs = 85.0;

	nstep = 5000;
	nplot = 333;	/* almost(!) every 50 seconds */

	printf("step  time   umin   umax     vmin   vmax     wmin   wmax     tmin   tmax      pmin   pmax\n");
	printf("-----------------------------------------------------------------------------------------------\n");

	/* set the initial conditions (at step n = 0) */
	ic(t1, t2, p1, p2, u1, v1, w1, u2, v2, w2, rho_t, rho_w, dx, dy, dz, I1, I2X, I2Y, I2Z, I2X_U, I2Y_V, I2Z_W, NX, NY, NZ, uPrt, tPrt, tBar, g, Cp);
	bc(t1, p1, u1, v1, w1, I1, I2X, I2Y, I2Z, I2X_U, I2Y_V, I2Z_W, NX, NY, NZ);

	/* find the minima and maxima for all fields (at step n = 0) */
	stats(t1, p1, u1, v1, w1, I1, I2X, I2Y, I2Z, I2X_U, I2Y_V, I2Z_W, NX, NY, NZ, 0, dt, tBar);
	
	/* plot the initial conditions (at step n = 0) */
	#pragma omp parallel for private(i, j, k) shared(t1, p1, u1, v1, w1, tPlot, pPlot, uPlot, vPlot, wPlot)
	for (i = I1; i <= I2X; i++)
		for (j = I1; j <= I2Y; j++)
			for (k = I1; k <= I2Z; k++)
			{
				tPlot[i - I1][j - I1][k - I1] = t1[i][j][k] - tBar;
				pPlot[i - I1][j - I1][k - I1] = p1[i][j][k];
				uPlot[i - I1][j - I1][k - I1] = (u1[i][j][k] + u1[i+1][j][k]) / 2.0;
				vPlot[i - I1][j - I1][k - I1] = (v1[i][j][k] + v1[i][j+1][k]) / 2.0;
				wPlot[i - I1][j - I1][k - I1] = (w1[i][j][k] + w1[i][j][k+1]) / 2.0;
			}

	putfield("T", 0.0, tPlot, NX, NY, NZ);
	putfield("P", 0.0, pPlot, NX, NY, NZ);
	putfield("U", 0.0, uPlot, NX, NY, NZ);
	putfield("V", 0.0, vPlot, NX, NY, NZ);
	putfield("W", 0.0, wPlot, NX, NY, NZ);

	/* time evolution */
	for (n = 1; n <= nstep; n++)
	{
		if (n == 1)
			tstep = dt;
		else
			tstep = 2.0 * dt;

		/* set boundary conditions */
		bc(t1, p1, u1, v1, w1, I1, I2X, I2Y, I2Z, I2X_U, I2Y_V, I2Z_W, NX, NY, NZ);
		bc(t2, p2, u2, v2, w2, I1, I2X, I2Y, I2Z, I2X_U, I2Y_V, I2Z_W, NX, NY, NZ);

		/* calculate the next step */
		advect(t2, t3, u1, v1, w1, u2, v2, w2, u3, v3, w3, dx, dy, dz, dt, tstep, I1, I2X, I2Y, I2Z, I2X_U, I2Y_V, I2Z_W, NX, NY, NZ);
		diffuse(t2, t3, u1, v1, w1, u3, v3, w3, K, Kt, dx, dy, dz, dt, tstep, I1, I2X, I2Y, I2Z, I2X_U, I2Y_V, I2Z_W, NX, NY, NZ);
		pgf(t2, p1, p3, u3, v3, w3, rho_t, rho_w, dx, dy, dz, tstep, I1, I2X, I2Y, I2Z, I2X_U, I2Y_V, I2Z_W, NX, NY, NZ, tBar, g, Cs);

		/* update the arrays at the end of each time step */
		update(t1, t2, t3, p1, p2, p3, u1, v1, w1, u2, v2, w2, u3, v3, w3, I1, I2X, I2Y, I2Z, I2X_U, I2Y_V, I2Z_W, NX, NY, NZ, n);

		/* find the minima and maxima for all fields */
		stats(t2, p2, u2, v2, w2, I1, I2X, I2Y, I2Z, I2X_U, I2Y_V, I2Z_W, NX, NY, NZ, n, dt, tBar);

		/* at each plot interval... */
		if (n % nplot == 0)
		{
			/* ...plot the fields */
			#pragma omp parallel for private(i, j, k) shared(t2, p2, u2, v2, w2, tPlot, pPlot, uPlot, vPlot, wPlot)
			for (i = I1; i <= I2X; i++)
				for (j = I1; j <= I2Y; j++)
					for (k = I1; k <= I2Z; k++)
					{
						tPlot[i - I1][j - I1][k - I1] = t2[i][j][k] - tBar;
						pPlot[i - I1][j - I1][k - I1] = p2[i][j][k];
						uPlot[i - I1][j - I1][k - I1] = (u2[i][j][k] + u2[i+1][j][k]) / 2.0;
						vPlot[i - I1][j - I1][k - I1] = (v2[i][j][k] + v2[i][j+1][k]) / 2.0;
						wPlot[i - I1][j - I1][k - I1] = (w2[i][j][k] + w2[i][j][k+1]) / 2.0;
					}

			putfield("T", dt * (float) n, tPlot, NX, NY, NZ);
			putfield("P", dt * (float) n, pPlot, NX, NY, NZ);
			putfield("U", dt * (float) n, uPlot, NX, NY, NZ);
			putfield("V", dt * (float) n, vPlot, NX, NY, NZ);
			putfield("W", dt * (float) n, wPlot, NX, NY, NZ);
		}
	}

	exit;
}
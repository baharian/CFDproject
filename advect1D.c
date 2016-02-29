/*
 *  ATMS 502 / CSE 566 -- Fall, 2012
 *  pgm5:  3D nonlinear quasi-compressible flow
 *  SOHEIL BAHARIAN
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "defs.h"

/* ======================= advect1D ====================================
 * ADVECT1D advects forward along only one dimension by one time step.
 *
 * Arguments:
 *	s1				real array		values at current step
 *	s2				real array		values at next step
 *	vel				real array		values of the velocity field
 *	dx				real			grid spacing
 *	dt				real			time step
 *	i1, i2			integers		indices bounding array data
 *	nx				integer			number of grid points
 *	advection_type	char			'L' is Lax-Wendroff, 'U' is Upstream,
 *									'T' is Takacs, 'C' is Crowley 6th-order,
 *									'P' is Piecewise-Linear
 */
void advect1D(s1, s2, vel, dx, dt, i1, i2, nx, advection_type)
int i1, i2, nx;
char advection_type;
float s1[NXDIM], s2[NXDIM], vel[NX+1], dx, dt;
{
	int i;
	float courant, r, ds, F[NX+1];

	switch (advection_type)
	{
		/* Lax-Wendroff method */
		case 'L':
			for (i = i1; i <= i2; i++)
			{
				courant = ((vel[i - i1] + vel[i - i1 + 1]) / 2.0) * (dt/dx);

				s2[i] = s1[i] - (courant / 2.0) * (s1[i+1] - s1[i-1]) + (powf(courant, 2.0) / 2.0) * (s1[i+1] - 2.0 * s1[i] + s1[i-1]);
			}
			
			break;


		/* Upsteam method */
		case 'U':
			for (i = i1; i <= i2; i++)
			{
				courant = ((vel[i - i1] + vel[i - i1 + 1]) / 2.0) * (dt/dx);

				if (courant > 0.0)
					s2[i] = s1[i] - courant * (s1[i] - s1[i-1]);
				else
					s2[i] = s1[i] - courant * (s1[i+1] - s1[i]);
			}

			break;


		/* Takacs method */
		case 'T':
			for (i = i1; i <= i2; i++)
			{
				courant = ((vel[i - i1] + vel[i - i1 + 1]) / 2.0) * (dt/dx);

				if (courant > 0.0)
					s2[i] = s1[i] - (courant / 2.0) * (s1[i+1] - s1[i-1]) + (powf(courant, 2.0) / 2.0) * (s1[i+1] - 2.0 * s1[i] + s1[i-1])
							- ((1.0 + fabs(courant)) / 6.0) * courant * (courant - 1.0) * (s1[i+1] - 3.0 * s1[i] + 3.0 * s1[i-1] - s1[i-2]);
				else
					s2[i] = s1[i] - (courant / 2.0) * (s1[i+1] - s1[i-1]) + (powf(courant, 2.0) / 2.0) * (s1[i+1] - 2.0 * s1[i] + s1[i-1])
							- ((1.0 + fabs(courant)) / 6.0) * courant * (courant + 1.0) * (s1[i-1] - 3.0 * s1[i] + 3.0 * s1[i+1] - s1[i+2]);
			}

			break;


		/* Crowley 6th-order method */
		case 'C':
			for (i = i1; i <= i2; i++)
			{
				courant = ((vel[i - i1] + vel[i - i1 + 1]) / 2.0) * (dt/dx);

				s2[i] = s1[i]
						+ (courant / 60.0) * (s1[i-3] - 9.0 * s1[i-2] + 45.0 * s1[i-1] - 45.0 * s1[i+1] + 9.0 * s1[i+2] - s1[i+3])
						+ (powf(courant, 2.0) / 360.0) * (2.0 * s1[i-3] - 27.0 * s1[i-2] + 270.0 * s1[i-1] - 490.0 * s1[i] + 270.0 * s1[i+1] - 27.0 * s1[i+2] + 2.0 * s1[i+3])
						+ (powf(courant, 3.0) / 48.0) * (- s1[i-3] + 8.0 * s1[i-2] - 13.0 * s1[i-1] + 13.0 * s1[i+1] - 8.0 * s1[i+2] + s1[i+3])
						+ (powf(courant, 4.0) / 144.0) * (- s1[i-3] + 12.0 * s1[i-2] - 39.0 * s1[i-1] + 56.0 * s1[i] - 39.0 * s1[i+1] + 12.0 * s1[i+2] - s1[i+3])
						+ (powf(courant, 5.0) / 240.0) * (s1[i-3] - 4.0 * s1[i-2] + 5.0 * s1[i-1] - 5.0 * s1[i+1] + 4.0 * s1[i+2] - s1[i+3])
						+ (powf(courant, 6.0) / 720.0) * (s1[i-3] - 6.0 * s1[i-2] + 15.0 * s1[i-1] - 20.0 * s1[i] + 15.0 * s1[i+1] - 6.0 * s1[i+2] + s1[i+3]);
			}

			break;


		/* Piecewise-Linear method */
		case 'P':
			#pragma omp parallel for private(i, r, ds) shared(s1, vel, F)
			for (i = 0; i < nx + 1; i++)
			{
				r = (dt/dx) * fabs(vel[i]);	/* r is calculated at the location of the staggered velocity */

				if (vel[i] < 0.0)
				{
					ds = (s1[i + i1 + 1] - s1[i + i1 - 1]) / 2.0;
					F[i] = r * (- s1[i + i1] + ((1.0 - r) / 2.0) * ds);
				}
				else
				{
					ds = (s1[i + i1] - s1[i + i1 - 2]) / 2.0;
					F[i] = r * (s1[i + i1 - 1] + ((1.0 - r) / 2.0) * ds);
				}
			}

			#pragma omp parallel for private(i) shared(s1, s2, vel, F)
			for (i = i1; i <= i2; i++)
				s2[i] = s1[i] - (F[i - i1 + 1] - F[i - i1]) + (dt/dx) * s1[i] * (vel[i - i1 + 1] - vel[i - i1]);

			break;


		default:
			printf("Integrate: Error! Unrecognized advection type: '%c'\n", advection_type);
			exit(1);
	}

	return;
}
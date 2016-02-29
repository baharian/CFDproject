/*
 *  ATMS 502 / CSE 566 -- Fall, 2012
 *  pgm5:  3D nonlinear quasi-compressible flow
 *  SOHEIL BAHARIAN
 */

/* definitions */
#define NX			300
#define NY			300
#define NZ			75

#define BC_WIDTH	2

#define I1			BC_WIDTH

#define I2X			I1 + NX - 1
#define I2Y			I1 + NY - 1
#define I2Z			I1 + NZ - 1

#define I2X_U		I1 + (NX + 1) - 1
#define I2Y_V		I1 + (NY + 1) - 1
#define I2Z_W		I1 + (NZ + 1) - 1

#define NXDIM		NX + 2 * BC_WIDTH
#define NYDIM		NY + 2 * BC_WIDTH
#define NZDIM		NZ + 2 * BC_WIDTH

#define NXDIM_U		(NX + 1) + 2 * BC_WIDTH
#define NYDIM_V		(NY + 1) + 2 * BC_WIDTH
#define NZDIM_W		(NZ + 1) + 2 * BC_WIDTH

#define	prt			2	/* number of perturbation */
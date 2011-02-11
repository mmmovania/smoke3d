/*
 *  smoke3D.cpp
 *  smoke3D
 *
 */

#include "types.h"
#include "utility.h"
#include "smoke3D.h"
#include "utility.h"
#include "solver.h"
#include "advect.h"
#include "render.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define	N			32
#define LIMIT		100
#define DT			0.1
#define SPHERE_R	0.2

FLOAT *** u[3] = { NULL, NULL, NULL };
FLOAT *** b = NULL;
FLOAT *** c = NULL;
int frame = 0;

void smoke3D::init() {
	// Allocate Memory
	u[0] = alloc3D(N+1,N,N);
	u[1] = alloc3D(N,N+1,N);
	u[2] = alloc3D(N,N,N+1);
	c = alloc3D(N,N,N);
	b = alloc3D(N,N,N);
	
	FOR_EVERY_X_FLOW {
		u[0][i][j][k] = 0.0;
	} END_FOR
	
	FOR_EVERY_Y_FLOW {
		u[1][i][j][k] = 0.0;
	} END_FOR
	
	FOR_EVERY_Z_FLOW {
		u[2][i][j][k] = 0.0;
	} END_FOR
	
	// Mark Wall Inside A Sphere
	int w = SPHERE_R*N;
	for( int i=-w; i<=w; i++ ) {
		for( int j=-w; j<=w; j++ ) {
			for( int k=-w; k<=w; k++ ) {
				if( hypot(hypot(i,j),k) < w ) {
					b[i+N/2][j+N/2][k+N/2] = 1.0;
				}
			}
		}
	}

#if _OPENMP
	printf( "OpenMP Detected.\n" );
#endif
}

static void enforce_boundary() {
	// Set Boundary Velocity Zero
	FOR_EVERY_X_FLOW {
		if( i==0 || i==N ) u[0][i][j][k] = 0.0;
	} END_FOR
	
	FOR_EVERY_Y_FLOW {
		if( j==0 || j==N ) u[1][i][j][k] = 0.0;
	} END_FOR
	
	FOR_EVERY_Z_FLOW {
		if( k==0 || k==N ) u[2][i][j][k] = 0.0;
	} END_FOR
	
	// Set Initial Smoke
	if( frame < LIMIT/2 ) {
		int w = N/7;
		for( int i=-w; i<=w; i++ ) {
			for( int j=-w; j<=w; j++ ) {
				if( hypot(i,j) < w ) {
					for( int k=0; k<6; k++ ) {
						u[1][(int)(N/2)+i][k+1][(int)(N/2)+j] = 0.0;
						c[(int)(N/2)+i][k+1][(int)(N/2)+j] = 1.0;
					}
				}
			}
		}
	}
	
	// Set Boundary Flow Around Sphere Zero
	FOR_EVERY_CELL {
		if( b[i][j][k] ) {
			c[i][j][k] = 0.0;
			u[0][i][j][k] = u[0][i+1][j][k] = 0.0;
			u[1][i][j][k] = u[1][i][j+1][k] = 0.0;
			u[2][i][j][k] = u[2][i][j][k+1] = 0.0;
		}
	} END_FOR
	
	// Give Some External Force
	FOR_EVERY_CELL {
		u[1][i][j][k] += 0.1*c[i][j][k];
	} END_FOR
}

static void project() {
	// Cell Width
	FLOAT h = 1.0/N;
	
	// Memory Allocation
	static FLOAT *** div = alloc3D(N,N,N);
	static FLOAT *** p = alloc3D(N,N,N);
	
	// Compute Divergence
	FOR_EVERY_CELL {
		div[i][j][k] = (u[0][i+1][j][k]-u[0][i][j][k]+
						u[1][i][j+1][k]-u[1][i][j][k]+
						u[2][i][j][k+1]-u[2][i][j][k]) / h;
	} END_FOR
	
	// Solve Pressure
	solver::solve( p, div, b, N );
	
	// Subtract Pressure Gradient
	FOR_EVERY_CELL {
		if( i>0 && i<N ) u[0][i][j][k] -= (p[i][j][k]-p[i-1][j][k])/h;
		if( j>0 && j<N ) u[1][i][j][k] -= (p[i][j][k]-p[i][j-1][k])/h;
		if( k>0 && k<N ) u[2][i][j][k] -= (p[i][j][k]-p[i][j][k-1])/h;
	} END_FOR
}

static void advection() {
	// Advect
	advect::advect( u, c, N, DT );
}

void smoke3D::simulateStep() {
	enforce_boundary();
	project();
	advection();
	render::render(c,SPHERE_R,N,frame);
	printf( "wrote frame %d\n", frame );
	frame ++;
	if( frame > LIMIT ) exit(0);
}

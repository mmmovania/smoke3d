/*
 *  render.cpp
 *  smoke3D
 *
 */

#include "render.h"
#include "utility.h"
#include "write_bmp.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define FOR_EVERY_PIXEL(W)	for( int ci=0; ci<W*W; ci++ ) { int i=ci%W; int j=ci/W;

static void normalize( FLOAT v[3] ) {
	FLOAT len = hypotf(v[0],hypotf(v[1],v[2]));
	if( len ) {
		v[0] /= len;
		v[1] /= len;
		v[2] /= len;
	}
}

static FLOAT interp( FLOAT *** d, int width, int height, int depth, FLOAT x, FLOAT y, FLOAT z )
{
	int i0, j0, k0, i1, j1, k1;
	FLOAT s0, t0, w0, s1, t1, w1;
	
	if (x<0.0) x=0.0; if (x>width) x=width; i0=min(width-2,(int)x); i1=i0+1;
	if (y<0.0) y=0.0; if (y>height) y=height; j0=min(height-2,(int)y); j1=j0+1;
	if (z<0.0) z=0.0; if (z>depth) z=depth; k0=min(depth-2,(int)z); k1=k0+1;
	
	s1 = x-i0; s0 = 1-s1; t1 = y-j0; t0 = 1-t1; w1 = z-k0; w0 = 1-w1;
	
	return		w0*(s0*(t0*d[i0][j0][k0]+t1*d[i0][j1][k0])+
				s1*(t0*d[i1][j0][k0]+t1*d[i1][j1][k0])) +
				w1*(s0*(t0*d[i0][j0][k1]+t1*d[i0][j1][k1])+
				s1*(t0*d[i1][j0][k1]+t1*d[i1][j1][k1]));
}

static double sample3D( FLOAT ***d, int N, FLOAT p[3] ) {
	if( p[0] < 0.0 || p[0] > 1.0 || 
		p[1] < 0.0 || p[1] > 1.0 ||
		p[2] < 0.0 || p[2] > 1.0 ) return 0.0;
	
	return interp( d, N, N, N, N*p[0], N*p[1], N*p[2] );
}

// Referenced site: https://mmack.wordpress.com/2010/11/01/adventures-in-fluid-simulation/
void render::render( FLOAT ***d, FLOAT sphere_r, int N, int frame ) {
	
	// Eye Position
	const FLOAT eyePos[3] = { 0.5, 0.5, -1.0 };
	
	// Light Position
	const FLOAT lightPos[3] = { 0.5, 1.5, 0.2 };
	
	// Light Intensity
	const FLOAT lightIntensity = 8.0;
	
	// Absortion
	const FLOAT absorption = 11.0;
	
	// Ray Samples
	const int numSamples = 128;
	
	// Ray Light Samples
	const int numLightSample = 64;
	
	// Maximum Distance
	const FLOAT maxDist = 3.0;
	
	// Ray Sampling Space
	const FLOAT stride = maxDist/numSamples;
	
	// Light Ray Sampling Space
	const FLOAT lstride = maxDist/numLightSample;
	
	// Window Size
	const int w = 256; 
	
	// Memory Allocation
	static unsigned char *image = new unsigned char[w*w*4];
	
	OPENMP_FOR FOR_EVERY_PIXEL(w) {
		// Transmittance
		FLOAT T = 1.0;
		
		// In-scattered Radiance
		FLOAT Lo = 0.0;
		
		// Pixel Position
		FLOAT pixPos[3] = {i/(FLOAT)w,j/(FLOAT)w};
		
		// Eye Vector
		FLOAT eyeVec[3] = {pixPos[0]-eyePos[0],pixPos[1]-eyePos[1],pixPos[2]-eyePos[2]};
		normalize(eyeVec);
		
		// Sphere Hit
		bool hitSphere = false;
		FLOAT dotProduct = 0.0;
		
		// Sample Density...
		for( int n=0; n<numSamples; n++ ) {
			
			// Sample Point
			FLOAT pos[3] = { eyePos[0]+stride*n*eyeVec[0], eyePos[1]+stride*n*eyeVec[1], eyePos[2]+stride*n*eyeVec[2] };
			
			// Sample Density
			FLOAT density = sample3D(d, N, pos);
			
			// Skip Empty Density
			if( density > 0.0 ) {
				
				// Attenuate
				T *= 1.0-density*stride*absorption;
				if( T <= 0.01 ) break;
				
				// Compute Light Ray
				FLOAT lightVec[3] = {lightPos[0]-pos[0],lightPos[1]-pos[1],lightPos[2]-pos[2]};
				normalize(lightVec);
				
				// Transmittance Along Light Ray
				FLOAT Tl = 1.0;
				
				// Sample Density Again...
				for( int m=1; m<numLightSample; m++ ) {
					
					// Sample Point
					FLOAT lpos[3] = { pos[0]+lstride*m*lightVec[0], pos[1]+lstride*m*lightVec[1], pos[2]+lstride*m*lightVec[2] };
					
					// Sphere Hit Test
					if( hypot(lpos[0]-0.5,hypot(lpos[1]-0.5,lpos[2]-0.5)) < sphere_r ) {
						Tl *= 1.0-exp(-3.0*lstride*m);
						break;
					}
					
					// Sample Density
					FLOAT ldensity = sample3D(d, N, lpos);
					
					// Attenuate
					Tl *= 1.0-absorption*lstride*ldensity;
					if( Tl <= 0.01 ) break;
				}
				
				FLOAT Li = lightIntensity*Tl;
				Lo += Li*T*density*stride;
			}
			
			// Sphere Collision ?
			if( hypot(pos[0]-0.5,hypot(pos[1]-0.5,pos[2]-0.5)) < sphere_r ) {
				
				// Compute Light Ray
				FLOAT lightVec[3] = {lightPos[0]-pos[0],lightPos[1]-pos[1],lightPos[2]-pos[2]};
				normalize(lightVec);
				
				// Compute Sphere Normal
				FLOAT normal[3] = { pos[0]-0.5, pos[1]-0.5, pos[2]-0.5 };
				normalize(normal);
				
				// Compute Dot Product
				dotProduct = max(0.1,normal[0]*lightVec[0] + normal[1]*lightVec[1] + normal[2]*lightVec[2]);
				hitSphere = true;
				break;
			}
		}
		
		// At Floor
		FLOAT Tf = 0.0;
		if( eyeVec[1] < 0.0 && hitSphere == false ) {
			FLOAT flen = -pixPos[1]/eyeVec[1];
			Tf = exp(-0.3*flen);
			
			// Compute Floor Intersection
			FLOAT pos[3] = {pixPos[0]+flen*eyeVec[0],pixPos[1]+flen*eyeVec[1],pixPos[2]+flen*eyeVec[2]};
			
			// Compute Light Ray
			FLOAT lightVec[3] = {lightPos[0]-pos[0],lightPos[1]-pos[1],lightPos[2]-pos[2]};
			normalize(lightVec);
			
			// Sample Density
			for( int m=1; m<numLightSample; m++ ) {
				
				// Sample Point
				FLOAT lpos[3] = { pos[0]+lstride*m*lightVec[0], pos[1]+lstride*m*lightVec[1], pos[2]+lstride*m*lightVec[2] };
				
				// Sample Density
				FLOAT ldensity = sample3D(d, N, lpos);
				
				// Attenuate
				Tf *= 1.0-0.5*absorption*lstride*ldensity;
				if( Tf <= 0.01 ) break;
			}
		}
		
		// Floor Color
		unsigned char floor_color[3] = { 75, 60, 45 };
		
		// Sphere Color
		unsigned char sphere_color[3] = { 50, 100, 150 };
		
		for( int k=0; k<3; k++ ) {
			image[4*(i+j*w)+k] = max(0,min(255,255*Lo + T*(Tf*floor_color[k] + dotProduct*sphere_color[k])));
		}
	} END_FOR
	
	char filename[64];
	sprintf( filename, "render_%d.bmp", frame );
	write_bmp( filename, image, w, w, false );
}
















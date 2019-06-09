#ifndef ANIM_UTILS
#define ANIM_UTILS

#include <stdio.h>

#define RADTODEG(radian) ((radian) * (180.0f/3.14159265358979323846264338327f))
#define DEGTORAD(degree) ((degree) * (3.14159265358979323846264338327f/180.0f))

// Convert HSL colour to RGB colour
void hslToRgb(double hRad, double s, double l, int & r, int & g, int &b);

// Convert RGB colour to HSL colour
void rgbToHsl(int r, int g, int b, double & h, double &s, double &l);

/**
 * Helper functions to print state to stream (typically for gnuplot)
 */
void fprintv
(
	FILE* stream,              // output stream
	// model parameters
	const   size_t N,            // number of particles
	// variables/buffers
	const   double* const  h,    // angles, in range [-pi,pi]
	const   double* const  x,    // x coords, in range [0,L]
	const   double* const  y     // y coords, in range [0,L]
);

void fprintv
(
	FILE* stream,              // output stream
								 // model parameters
	const   size_t N,            // number of particles
								 // variables/buffers
	const   double* const  h,    // angles, in range [-pi,pi]
	const   double* const  x,    // x coords, in range [0,L]
	const   double* const  y,     // y coords, in range [0,L]
	const	int* const r,
	const	int* const g,
	const	int* const b
);

void fprintv
(
	FILE* stream,              // output stream
								 // model parameters
	const   size_t N,            // number of particles
								 // variables/buffers
	const   double* const  x,    // x coords, in range [0,L]
	const   double* const  y     // y coords, in range [0,L]
);

#endif
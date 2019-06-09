#include "../include/anim_utils.hpp"

#include "../include/defs.hpp"
#include <fmt/format.h>
#include <algorithm>


double hue2rgb(double p, double q, double t)
{
	if(t < 0) t += 1;
	if(t > 1) t -= 1;
	if(t < 1.0f / 6) return p + (q - p) * 6 * t;
	if(t < 1.0f / 2) return q;
	if(t < 2.0f / 3) return p + (q - p) * (2.0f / 3 - t) * 6;
	return p;
}

//hRad = [0,1]
void hslToRgb(double hRad, double s, double l, int & r, int & g, int &b)
{
	double rd, gd, bd;
	double h = RADTODEG(hRad * TWOPI) / 360.0;

	if(s == 0) {
		rd = gd = bd = l; // achromatic
	}
	else {
		double q = l < 0.5 ? l * (1 + s) : l + s - l * s;
		double p = 2 * l - q;
		rd = hue2rgb(p, q, h + 1.0f / 3);
		gd = hue2rgb(p, q, h);
		bd = hue2rgb(p, q, h - 1.0f / 3);
	}

	r = (int)(rd * 255);
	g = (int)(gd * 255);
	b = (int)(bd * 255);
}

void rgbToHsl(int r, int g, int b, double & h, double &s, double &l)
{
	float rf = (float)r / 255.0f;
	float gf = (float)g / 255.0f;
	float bf = (float)b / 255.0f;
	float max = std::max(rf, std::max(gf, bf));
	float min = std::min(rf, std::min(gf, bf));
	h = (max + min) * 0.5f;
	s = (max + min) * 0.5f;
	l = (max + min) * 0.5f;

	if(max == min)
	{
		h = s = 0; // achromatic
	}
	else
	{
		float d = max - min;
		s = l > 0.5 ? d / (2 - max - min) : d / (max + min);
		if(max == rf)
			h = (gf - bf) / d + (gf < bf ? 6 : 0);
		else if(max == gf)
			h = (bf - rf) / d + 2;
		else //if(max == rgb.b)
			h = (rf - gf) / d + 4;
		h /= 6;
	}
}

void fprintv
(
	FILE* stream,              // output stream
	// model parameters
	const   size_t N,            // number of particles
	// variables/buffers
	const   double* const  h,    // angles, in range [-pi,pi]
	const   double* const  x,    // x coords, in range [0,L]
	const   double* const  y     // y coords, in range [0,L]
)
{
	for (size_t i = 0; i < N; ++i) { // for each particle
		fmt::print(stream, "{:<20.16} {:<20.16} {:<20.16}\n", h[i], x[i], y[i]);
	}
}

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

)
{
	for (size_t i = 0; i < N; ++i) { // for each particle
		fmt::print(stream, "{:<20.16} {:<20.16} {:<20.16} {} {} {}\n", h[i], x[i], y[i], r[i], g[i], b[i]);
	}
}

void fprintv
(
	FILE* stream,              // output stream
								 // model parameters
	const   size_t N,            // number of particles
								 // variables/buffers
	const   double* const  x,    // x coords, in range [0,L]
	const   double* const  y     // y coords, in range [0,L]
)
{
	for (size_t i = 0; i < N; ++i) { // for each particle
		fmt::print(stream, "{:<20.16} {:<20.16}\n", x[i], y[i]);
	}
}
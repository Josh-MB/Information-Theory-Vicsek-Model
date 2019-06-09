#ifndef UTILS_H
#define UTILS_H

// macros and low-level inline functions

#include <clara.hpp>
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <vector>
#include <queue>
#include <functional>
#include <random>

// test for IEEE 754 compliance (let's you use memset to zero out fp arrays and much more)
#ifdef ___LINUX___
#ifndef __STDC_IEC_559__
#error "Hmm, not IEEE 754 compliant ... trouble ahead"
#endif
#endif

enum SIM
{
	MIBIN = 0,
	TEBIN,
	GTEBIN,
	PARAMS, //Used for writing out order param/suscept etc when no sims are performed
	SIM_COUNT,
};

// silly, but correct :)
#define PI     (double)3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348L
#define TWOPI  (double)6.2831853071795864769252867665590057683943387987502116419498891846156328125724179972560696L
#define ITWOPI (double)0.1591549430918953357688837633725143620344596457404564487476673440588967976342265350901138L

static const size_t GK = 6; // Number of Y dimensions to use for GTE

//On-lattice wrapping
inline size_t wrap_minus(size_t u, size_t L1) { return (u == 0) ? L1 : u - 1; }; // -1 with lattice wrap
inline size_t wrap_plus(size_t u, size_t L1) { return (u == L1) ? 0 : u + 1; }; // +1 with lattice wrap

inline double cwrap(const double x, const double L) // return coord wrapped to [0,L)
{
	return x-L*floor(x/L);
}

inline double rwrap(const double x, const double L) // return coord wrapped to range [-L/2, L/2)
{
	return x - L*floor((x / L) + 0.5);
}

inline double dwrap(const double d, const double L) // return distance in [0,L) wrapped to [0,L/2]
{
	const double dw = fabs(d);
	return dw < L/2 ? dw : L-dw;
	//Almost same as above line, except returns -dw when dw < L/2. Which is not a problem since it is only used in dsqwrap and it's squared, so sign doesn't matter
	//return (dw >= L*0.5)*L-dw;
	//return (dw >= L*0.5)*(L-2*dw)+dw;
}

//Returns the offset from x0 to x1, such that cwrap(x0+offset, L) == x1
inline double offsetwrap(const double x0, const double x1, const double L)
{
	const double dw = fabs(x1 - x0);
	return dw < L / 2 ? x1 - x0 : (x0 >= x1 ? L - dw : dw - L);
}

inline double dsqwrap(const double dx, const double dy, const double L) // return planar squared distance in [0,L) x [0,L) wrapped to [0,L/2] x [0,L/2]
{
	const double dxw = dwrap(dx,L);
	const double dyw = dwrap(dy,L);
	return dxw*dxw+dyw*dyw;
}

inline double awrap(const double h) // return angle wrapped to [-pi,pi)
{
	return h-TWOPI*floor(ITWOPI*h+0.5);
}

inline size_t binit(const double h, const size_t B) // return bin number of angle in range 0 .. B-1. IMPORTANT: ASSUMES ANGLE WRAPPED to [-pi,pi)
{
	const long b = (long)(floor((double)B*(ITWOPI*h+0.5)));
	return (size_t)(b < 0 ? 0 : b >= (long)B ? (long)B-1 : b); // note: could have b < 0 if h = -pi - eps, or b >= B if h =  pi + eps
	//return (size_t)((b >= 0) * ((b >= (long)B) * ((long)B - 1) + (!(b >= (long)B)) * b));//No faster really
	//A fraction faster than above two implementations. But assumes that
	//B >= 0 in all cases. Need to confirm.
	//return (size_t)((b >= (long)B)*((long)B-1) + (!(b >= (long)B))*(b>=0)*b);
}

inline size_t binit2(const double h, const size_t B) // return bin number of angle in range 0 .. B-1. IMPORTANT: ASSUMES ANGLE WRAPPED to [-2pi,2pi)
{
	const long b = (long)(floor((double)B*((0.5*ITWOPI)*h+0.5)));
	return (size_t)(b < 0 ? 0 : b >= (long)B ? (long)B-1 : b);
}

inline size_t binit3(const double h, const size_t B) // return bin number of angle in range 0 .. B-1. IMPORTANT: ASSUMES ANGLE WRAPPED to [0,pi)
{
	const long b = (long)(floor((double)B*((2*ITWOPI)*h)));
	return (size_t)(b < 0 ? 0 : b >= (long)B ? (long)B-1 : b);
}

inline double mxlogx(const double x) // -x ln(x), or 0 if x = 0
{
	return x < -DBL_EPSILON ? NAN : x <  DBL_EPSILON ? 0.0 : -x*log(x);
}

//Calc cross product assuming u_z = v_z = 0
inline double crossProduct_vxv(const double ux, const double uy, const double vx, const double vy)
{
	return ux*vy - uy*vx;
}

//Calc cross product assuming u = (0,0,u_z) and v=(v_x, v_y, 0), where
//s=uxv, and s = (s_x, s_y, 0)
inline void crossProduct_sxv(const double uz, const double vx, const double vy, double& sx, double& sy)
{
	sx = -uz*vy;
	sy =  uz*vx;
}

typedef std::pair<int, int> InteractionPair;

struct ksg_point {
	double h1, h2, h1new, h2Agg;
};


// Draws (without replacement) M elements from indices to shuffle.
// These M elements are placed at the front of indicesToShuffle (all elements
// remain in array, just now changed order).
// Will touch maximum 2*M points, much more efficient than using a complete random_shuffle
// which would touch all points.
// Note, there is a slight bias for lower numbers, however, since we are using a 64 bit
// mersenne twister and generating numbers betweon 0 and 2^64-1, this bias just means the first
// indices.size() % 2^64 elements will appear indices.size() / 2^64 percent more often than later
// numbers. For indices.size of ~268mil (2^28), this is only a 1x10^-11% increase.
inline void shuffleFirstMIndices(std::vector<InteractionPair>& indicesToShuffle,
						  size_t M,		// Randomise first M indices out of whole set
						  std::mt19937_64& rng            // 64-bit Mersenne Twister
)
{
	size_t left = indicesToShuffle.size();
	size_t begin = 0;
	while(M--) {
		size_t r = begin;
		r += rng() % left;
		std::swap(indicesToShuffle[begin], indicesToShuffle[r]);
		++begin;
		--left;
	}
}

// Calculates the consensus vector for each particle
// while only using the closest n neighbours rather than
// all those within interaction radius
void restrictedConsensusVector(const size_t N,
							 const double L,
							 const double* const x,
							 const double* const y,
							 const double* const h,
							 double* const hAgg,
							 const size_t numNeighbours);

// Parameters used during the initialisation of model
struct GeneratorParams {
	uint64_t seed{ 0 };
	size_t N{ 500 };
	size_t S{ 1000 };
	size_t topo_neighbours{ 6 };
	double rho{ 0.25 };
	double v{ 0.3 };
	double eta{ 0.3 };
	double L{ sqrt(static_cast<double>(N) / rho) };
	double chi{ 1.25 };
	double J{ 0.8 };
	double viscosity{ 0.3 };
	double dt_factor{ 1.0 };
	double force_dt{ 0. };
	int imethod{ 0 };
	int umethod{ 0 };
	int discretise{ 0 };
	bool rotate_frame{ false };
};

// Parameters used for simulation and calculation of
// info theory metrics
struct SimParams {
	uint64_t seed{ 0 };
	size_t U{ 1000 };
	size_t B{ 16 };
	size_t GKavg{ 0 };
	size_t KSG_neighbours{ 3 };
	double subset_size{ 0 };
	int te_shuffle_dim{ 1 };
	int subset_runs{ 1 };
	int ksg_local{ 0 };
	int record_T_steps{ 1 };
	int ksg_gte_dims{ 2 }; // Number of dimensions to use for KSG calculation of GTE. 0 = GK+2, 1 = unused, 2 = Tgl^2D, 3 = Tgl^3D
	int hist_gte_dims{ 1 }; // Number of dimensions to use for histogram estimation of GTE, 0=disabled, or 1 or 2.
	bool sims[SIM_COUNT]{ false,false,false,false };
	bool angleWrite{ false };
	bool reseed{ false };
	bool ksg_t{ false };
	bool ksg_w{ false };
	bool shuffle{ false };
};

// Command line argument parsers
clara::detail::Parser generatorParamsParser(GeneratorParams& gp);
clara::detail::Parser simParamsParser(SimParams & sp);
#endif // UTILS_H

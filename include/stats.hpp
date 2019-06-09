#ifndef STATS_H
#define STATS_H

#include "defs.hpp"
#include <map>
#include <vector>
#include <random>

void order_param(
	const   size_t        N,     // number of particles
	const   double* const h,     // angles, in range [-pi,pi]
	double* const         zx,    // order parmeter x coord
	double* const         zy,     // order parmeter y coord
	double  const		  cumulRot = 0 //Cumulative frame rotation
);

void nn_histogram(
	const   size_t        N,     // number of particles
	const   double        L,     // linear size
	const   double* const x,     // x coords, in range [0,L]
	const   double* const y,     // y coords, in range [0,L]
	uint64_t*  const         count  // number of neighbours histogram (accumulate)
);

void MI_histogram(
	const  size_t        N,      // number of particles
	const  double* const h,      // angles, in range [-pi,pi]
	const  double        cumulativeRotation,
	const  size_t        B,      // number of bins
	uint64_t* const         count,  // histogram (accumulate)
	size_t*  const	     bins,   // bins of every particle 
	std::vector<int> const & interactionSet
);

void MI_1D_histogram(
	const  size_t        N,      // number of particles
	const  double        L,      // linear size
	const  double* const h,      // angles, in range [-pi,pi]
	const  double* const x,      // x coords, in range [0,L]
	const  double* const y,      // y coords, in range [0,L]
	const  size_t        B,      // number of bins
	uint64_t* const         count,  // 1-dim histogram (accumulate)
	uint64_t* const		 marginalCount
);

void GTE_1D_histogram(
	const  size_t        N,      // number of particles
	const  double* const h,      // angles, in range [-pi,pi]
	const  double* const hnew,
	const  size_t        B,      // number of bins
	uint64_t* const         count,  // 1-dim histogram (accumulate)
	const  double		 rotation //Rotation of last frame, such that \Theta' is in same reference frame as \Theta
);

void GTE_2D_histogram(
	const  size_t        N,      // number of particles
	const  double* const h,      // angles, in range [-pi,pi]
	const  double* const hnew,
	const  size_t        B,      // number of bins
	uint64_t* const         count,  // 2-dim histogram (accumulate)
	const  double		 rotation //Rotation of last frame, such that \Theta' is in same reference frame as \Theta
);

void MI_histogram_shuffled(
	const size_t 		N,
	const size_t		U,      // number of updates
	const double		L,
	const double* const 	h,
	const double* const	x,
	const double* const 	y,
	const size_t		B,
	uint64_t* const		count,
	size_t* const		bins,
	std::vector<std::vector<int> > const & fullInteractionSet,
	std::mt19937_64 &rng
);

void MI_histogram_all_pairs(
	const  size_t        N,      // number of particles
	const  double* const h,      // angles, in range [-pi,pi]
	const  size_t        B,      // number of bins
	uint64_t* const         count,  // histogram (accumulate)
	size_t*  const	     bins    // bins of every particle 
);

void TE_histogram(
	const  size_t        N,      // number of particles
	const  double* const h,      // angles, in range [-pi,pi]
	const  double* const hnew,   // new angles
	const  double        rotation,
	const  double        cumulativeRotation,
	const  size_t        B,      // number of bins
	uint64_t* const         count,  // histogram (accumulate)
	size_t*  const	     bins,    // bins of every particle 
	std::vector<int> const & interactionSet
);

void TE_histogram_shuffled(
	const  size_t        N,      // number of particles
	const  size_t        U,      // number of updates
	const  double* const h,      // angles, in range [-pi,pi]
	const  size_t        B,      // number of bins
	uint64_t* const         count,  // histogram (accumulate)
	size_t*  const	     bins,   // bins of every particle 
	const  int	     d,	     // dimension to shuffle (0 = x, 1 = y, 2 = w)
	std::vector<std::vector<int> > const & fullInteractionSet,
	std::mt19937_64& rng
);

void insort(
	const  size_t   N,            // number of elements to sort
	double*         const dsq,    // elements to sort
	size_t*         const idx     // indices of elements
);

void GTE_histogram(
	const   size_t        N,      // number of particles
	const   double        L,      // linear size
	const   double* const h,      // angles, in range [-pi,pi]
	const   double* const hnew,   // new angles
	const   double* const x,      // x coords, in range [0,L]
	const   double* const y,      // y coords, in range [0,L]
	const   size_t        K,      // number of nearest neighbours in source
	double*         const dsq,    // square distances buffer
	size_t*         const idx,    // indices buffer
	const   size_t        B,      // number of bins
	uint64_t*  const         count,  // histogram (accumulate)
	size_t*  const	     bins    // bins of every particle 
);

void GTE_avg_histogram(
	const  size_t        N,      // number of particles
	const  double* const h,      // angles, in range [-pi,pi]
	const  double* const hnew,   // new angles
	const  double 	     rotation,
	const  double        cumulativeRotation,
	const  double* const hAgg,   // neighbour angles
	const  size_t        B,      // number of bins
	uint64_t* const         count  // histogram (accumulate)
);

void GTE_avg_histogram_shuffled(
	const  size_t        N,      // number of particles
	const  size_t        U,      // number of updates
	const  double* const h,      // angles, in range [-pi,pi]
	const  double* const hAgg,   // neighbour angles
	const  size_t        B,      // number of bins
	uint64_t* const         count,  // 3-dim histogram (accumulate)
	const  int	     d,	     // dimension to shuffle (0 = x, 1 = y, 2 = w)
	std::mt19937_64& rng
);

//Calculates the center of mass of the entire flock, accounting for periodic
//boundary conditions
void calcCenterMass(const size_t N, double const * const x, double const * const y, const double L,
					double &cx, double &cy);

void spatioTemporalCorrelation(const size_t N,
							   const double L,
							   double const * const x0,
							   double const * const y0,
							   double const * const h0,
							   double const * const xt,
							   double const * const yt,
							   double const * const ht,
							   double const cm_x0,
							   double const cm_y0,
							   double const cm_xt,
							   double const cm_yt,
							   double const rho,
							   double const v,
							   const size_t numR,
							   const double dr,
							   double const * const rVals,
							   double * const Crt_accum
);

void fourierSpatioTemporalCorrelation(const size_t N,
									  const double L,
									  double const * const x0,
									  double const * const y0,
									  double const * const h0,
									  double const * const xt,
									  double const * const yt,
									  double const * const ht,
									  double const cm_x0,
									  double const cm_y0,
									  double const cm_xt,
									  double const cm_yt,
									  double const v,
									  const size_t numK,
									  double const * const kVals,
									  double * const Crt_accum
);

class ITStatTracker {
private:
	double total = 0.;
	double min_ = 0.;
	double max_ = 0.;
	size_t count = 0;
public:
	void accumulate(double value) {
		total += value;
		min_ = std::min(min_, value);
		max_ = std::max(max_, value);
		++count;
	}
	double min() { return min_; }
	double max() { return max_; }
	double avg() { return total / static_cast<double>(count); }
};

struct AvgMinMax {
	double avg;
	double min;
	double max;
};
struct ITStats {
	ITStatTracker mi;
	ITStatTracker te;
	ITStatTracker gte;

	AvgMinMax operator[] (int pos) {
		switch (pos) {
		case MIBIN:
			return { mi.avg(), mi.min(), mi.max() };
		case TEBIN:
			return { te.avg(), te.min(), te.max() };
		case GTEBIN:
			return { gte.avg(), gte.min(), gte.max() };
		default:
			return { 0., 0., 0. };
		}
	}
};
#endif // STATS_H

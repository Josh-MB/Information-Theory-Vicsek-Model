#ifndef MODEL_H
#define MODEL_H

#include "../include/vpUtil.hpp"
#include "../include/vpTree.hpp"
#include "connectedFlocks.hpp"
#include "utils.hpp"

#include "defs.hpp"
#include <random>
#include <vector>

#define DECIMATE_INLINE

// Initialise flock state
void initialise(
	const   int    imode,        // angle initialisation mode
	// model parameters
	const   size_t N,            // number of particles
	const   double L,            // linear size
	// variables/buffers
	double* const  h,            // angles, in range [-pi,pi]
	double* const  x,            // x coords, in range [0,L]
	double* const  y,            // y coords, in range [0,L]
	double* const  s,            // y coords, in range [0,L]
	// random number generation
	std::mt19937_64& rng,           // 64-bit Mersenne Twister
	bool onLattice = false		 // Set to true to align particles to lattice layout (for XY model)
);

// Original vicsek model. Not used
void vicsek(
	// model parameters
	const   size_t N,            // number of particles
	const   double L,            // linear size
	const   double v,            // particle velocity
	const   double eta,          // noise parameter (stdandard deviation)
	// variables/buffers
	const   double* const  h,    // angles, in range [-pi,pi]
	const   double* const  x,    // x coords, in range [0,L]
	const   double* const  y,    // y coords, in range [0,L]
	double* const          hnew, // new angles
	double* const          xnew, // new x coords
	double* const          ynew, // new y coords
	// random number generation
	std::mt19937_64& rng,            // 64-bit Mersenne Twister
	const int updateMethod = 0    //0 - metric, 1 - topological
);

// Vicsek model using metric interactions
void vicsekMetric(
	// model parameters
	const   size_t N,            // number of particles
	const   double L,            // linear size
	const   double v,            // particle velocity
	const   double eta,          // noise parameter (stdandard deviation)
	// variables/buffers
	FlockState& fs,
	double* const		   hAgg, // Average of neighbours for each particle
	// random number generation
	std::mt19937_64& rng,            // 64-bit Mersenne Twister
	const   int updateMethod,  //0 = backwards, 1 = forwards
	VpTreeEuclidean2D & tree,
	std::vector<int>* interactionSet = 0,
	unsigned long int* interactionCount = 0,
	int discretise = 0,
	ConnectedSets* connectedSets = 0,
	bool rotate = false,
	double *rotation = nullptr
	);

// Vicsek metric interactions using multithreading
void vicsekMetricThreaded(
	// model parameters
	const   size_t N,            // number of particles
	const   double L,            // linear size
	const   double v,            // particle velocity
	const   double eta,          // noise parameter (stdandard deviation)
	// variables/buffers
	FlockState& fs,
	double* const		   hAgg, // Average of neighbours for each particle
	// random number generation
	std::mt19937_64& rng,            // 64-bit Mersenne Twister
	const   int updateMethod,  //0 = backwards, 1 = forwards
	VpTreeEuclidean2D & tree,
	std::vector<int>* interactionSet = 0,
	const int numThreads = 1
);

// 2D Inertial Spin Model. Differs from the full ISM in that particle
// headings are restricted to the XY plane, meaning spins are simply 
// a scalar (i.e. vector of the form (0, 0, z)) ensuring s.v = 0
void inertialSpinModel(
	// model parameters
	const   size_t N,            // number of particles
	const   double L,            // linear size
	const   double v,            // particle velocity
	const   double sigma,          // noise parameter (stdandard deviation)
	const   double chi,			 // Generalised moment of inertia (=1.25)
	const   double vis,			 // viscosity constant (\eta in Cavagna15, renamed. (=0.3)
	const   double J,			 // alignment strength (=0.8)
	const   double dt,			 // time delta
	const   size_t NC,			 // num neighbours
								 // variables/buffers
	FlockState& fs,				 // flock state. Angles, x, y. Angles in range [-pi,pi]
	const   double* const  s,	 // spin of particles
	double* const		   snew, // new spins
	double* const		   hAgg, // Average of neighbours for each particle
								 // random number generation
	std::mt19937_64& engine,	//Generator for gaussian rng
	std::normal_distribution<>& gaussian_rng, //Gaussian rng for ISM noise
	const   int updateMethod,  //0 = backwards, 1 = forwards
	VpTreeEuclidean2D & tree,
	double & hamiltonianEnergy,
	std::vector<int>* interactionSet = 0,
	unsigned long int* interactionCount = 0,
	ConnectedSets* connectedSets = 0
);

void vicsekTopological(
	// model parameters
	const   size_t N,            // number of particles
	const   double L,            // linear size
	const   double v,            // particle velocity
	const   double eta,          // noise parameter (stdandard deviation)
	const	double dt,			 // timestep
	// variables/buffers
	FlockState& fs,
	double* const	       hAgg, // Neighbours averaged into one
	// random number generation
	std::mt19937_64& rng,            // 64-bit Mersenne Twister
	const	size_t NC,			 //Number of topological neighbours
	VpTreeEuclidean2D & tree,
	const   int updateMethod,  //0 = backwards, 1 = forwards
	std::vector<int>* interactionSet = 0,
	unsigned long int* interactionCount = 0
	);

void xyModel(
	const size_t N,
	const double L,
	const double eta,
	FlockState& fs,
	double* const hAgg,
	std::mt19937_64& rng
);

void xy(
	// model parameters
	const   size_t N,            // number of particles
	const   double L,            // linear size
	const   double eta,          // noise parameter (stdandard deviation)
	// variables/buffers
	const   double* const  h,    // angles, in range [-pi,pi]
	const   double* const  x,    // x coords, in range [0,L]
	const   double* const  y,    // y coords, in range [0,L]
	double* const          hnew, // new angles
	double* const          xnew, // new x coords
	double* const          ynew, // new y coords
	// random number generation
	std::mt19937_64& rng            // 64-bit Mersenne Twister
);

void kuravic(
	// model parameters
	const   size_t N,            // number of particles
	const   double L,            // linear size
	const   double v,            // particle velocity
	const   double w,            // particle angular velocity
	const   double eta,          // noise parameter (stdandard deviation)
	// variables/buffers
	const   double* const  h,    // angles, in range [-pi,pi]
	const   double* const  x,    // x coords, in range [0,L]
	const   double* const  y,    // y coords, in range [0,L]
	double* const          hnew, // new angles
	double* const          xnew, // new x coords
	double* const          ynew, // new y coords
	// random number generation
	std::mt19937_64& rng           // 64-bit Mersenne Twister

);

// Rotate model frame
void rotateModel(
	const size_t N,
	const double L,
	const double rotation,
	const double* const h,
	const double* const x,
	const double* const y,
	double* const hnew,
	double* const xnew,
	double* const ynew
);

inline double jitter(const double vx, const double vy, const double eta, 
	std::mt19937_64& rng, std::uniform_real_distribution<double>& dis)
{
	return awrap(atan2(vy,vx) + eta*(dis(rng)-0.5));
}

#endif // MODEL_H

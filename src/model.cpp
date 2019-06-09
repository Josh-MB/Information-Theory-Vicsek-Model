#include "../include/model.hpp"
#include "../include/stats.hpp"

#include <fmt/format.h>
#include <stdlib.h>
#include <algorithm>
#include <string.h>
#include <queue>
#include <random>
#include <omp.h>

#ifdef ___OSX___
	#define sincos __sincos
#endif



void initialise(
	const   int    imode,        // angle initialisation mode
	// model parameters
	const   size_t N,            // number of particles
	const   double L,            // linear size
	// variables/buffers
	double* const  h,            // angles, in range [-pi,pi]
	double* const  x,            // x coords, in range [0,L]
	double* const  y,            // y coords, in range [0,L]
	double* const  s,            // spin, in range [-1,1]
	// random number generation
	std::mt19937_64&  rng,           // 64-bit Mersenne Twister
	bool onLattice		 // Set to true to align particles to lattice layout (for XY model)
)
{
	std::uniform_real_distribution<double> dis;
	if      (imode == 0) { // set all angle independently random
		for (size_t i=0;i<N;++i) h[i] = PI*(2.0*dis(rng)-1.0); // in range [-pi,pi]
	}
	else if (imode == 1) { // set all angles to some particular random angle
		const double hh = PI*(2.0*dis(rng)-1.0); // in range [-pi,pi]
		for (size_t i=0;i<N;++i) h[i] = hh;
	}
	else if (imode == 2) { // set all angles to zero
		for (size_t i=0;i<N;++i) h[i] = 0.0;
	}
	else if(imode == 3 || imode == 4) { // set all angles to zero
		for(size_t i = 0; i < N; ++i) h[i] = 0.0;
	}
	else EEXIT("unknown initialisation mode");
	//h[0] = 1.0;
	if (onLattice) {
		// Aligns particle position such that particle
		// N is located at (x,y) = (N % L, N / L)
		// Technically not actually needed, except for visualisation
		size_t L1 = static_cast<size_t>(L);
		for (size_t j = 0; j < L; ++j) {
			for (size_t i = 0; i < L; ++i) {
				auto idx = j * L1 + i;
				x[idx] = static_cast<double>(i);
				y[idx] = static_cast<double>(j);
			}
		}
	}
	else {
		if(imode == 3) {
			for(size_t i = 0; i < N; ++i) { // for each particle
				x[i] = 10* dis(rng)+15;
				y[i] = 10* dis(rng)+10;
				s[i] = 0.1*dis(rng) - 0.05;
				//s[i] = 0;
				//s[i] = 2*mt64_rand(rng) - 1;
			}
			x[0] = 20;
			y[0] = 15;
		}
		else if(imode == 4) {
			for(size_t i = 0; i < N; ++i) {
				double phi = dis(rng) * TWOPI;
				double u = dis(rng);
				double r = L * 0.5 * std::sqrt(u);

				x[i] = r * cos(phi) + L * 0.5;
				y[i] = r * sin(phi) + L * 0.5;
				s[i] = 0;
			}
		}
		else {
			for(size_t i = 0; i < N; ++i) { // for each particle
				x[i] = L* dis(rng);
				y[i] = L* dis(rng);
				//s[i] = 0.1*mt64_rand(rng) - 0.05;
				//s[i] = 0;
				s[i] = 0.1*dis(rng) - 0.05;
				//s[i] = 2 * mt64_rand(rng) - 1;
			}
		}
	}
}

void vicsekMetric(
	// model parameters
	const   size_t N,            // number of particles
	const   double L,            // linear size
	const   double v,            // particle velocity
	const   double eta,          // noise parameter (stdandard deviation)
	// variables/buffers
	FlockState &fs,
	double* const		   hAgg, // Average of neighbours for each particle
	// random number generation
	std::mt19937_64& rng,            // 64-bit Mersenne Twister
	const   int updateMethod,  //0 = backwards, 1 = forwards
	VpTreeEuclidean2D & tree,
	std::vector<int>* interactionSet,
	unsigned long int* interactionCount,
	int discretise,
	ConnectedSets* connectedSets,
	bool rotate,
	double *rotation
	)
{
	std::uniform_real_distribution<double> dis;
	auto h = fs.h.reader(), x = fs.x.reader(), y = fs.y.reader();
	auto hnew = fs.h.writer(), xnew = fs.x.writer(), ynew = fs.y.writer();

	std::vector<VPoint<2> > points(N);
	// precalculate unit vector velocities for every particle
	std::vector<double> xangles(N);
	std::vector<double> yangles(N);

	if(rotate) {
		if(!rotation) PEEXIT("Rotation parameter is nullptr");
		*rotation = TWOPI * (dis(rng));
		std::vector<double> rotX(N), rotY(N), hRot(N);
		rotateModel(N, L, *rotation, h, x, y, hRot.data(), rotX.data(), rotY.data());
		for(size_t i = 0; i < N; ++i) {
			sincos(hRot[i], &yangles[i], &xangles[i]);
			points[i].values[0] = rotX[i];
			points[i].values[1] = rotY[i];
			points[i].idx = (int)i;
		}
	}
	else {
		for(size_t i = 0; i < N; ++i) {
			sincos(h[i], &yangles[i], &xangles[i]);
			points[i].values[0] = x[i];
			points[i].values[1] = y[i];
			points[i].idx = (int)i;
		}
	}
	tree.create(points);

	std::vector<VPoint<2> > results;
	results.reserve(50);
	std::vector<int> indices;
	std::vector<double> distances;

	for(size_t i = 0; i<N; ++i) { // for each particle
		int id_i = points[i].idx;
		results.clear();
		indices.clear();
		distances.clear();
		//Work out neighbours and dists. Dnc = Distance to NCth neighbour
		int num = tree.fr_search(points[i], 1.0, &results);

		double xanv = 0.0, yanv = 0.0;
		for(int j = 0; j < num; ++j) {
			int id_j = results[j].idx;
			xanv += xangles[id_j];
			yanv += yangles[id_j];

			//add to interaction set
			if(interactionSet != 0 && id_i != id_j) {
				interactionSet->push_back(id_i);
				interactionSet->push_back(id_j);
			}
			if(interactionCount != 0 && id_i != id_j)
				(*interactionCount) += 2;

			if(connectedSets != 0 && id_i != id_j) {
				connectedSets->connect(id_i, id_j);
			}
		}
		if(connectedSets != 0 && num <= 1)
			connectedSets->add(id_i);

		// Store the aggregate of all interacting particles, without noise, for use with GTE
		hAgg[id_i] = awrap(atan2(yanv, xanv));
		hnew[id_i] = awrap(jitter(xanv, yanv, eta, rng, dis));

		if(discretise != 0) {
			double const d = discretise * ITWOPI;
			for(int idx = 0; idx < (int)N; ++idx) {
				hnew[idx] = std::floor(hnew[idx] * d) / (double)d;
			}
		}

		double xu, yu;
		if(updateMethod == 1) //Forwards
			sincos(hnew[id_i], &yu, &xu); // u is unit vector in direction hnew[i]
		else //Backwards
			sincos(h[id_i], &yu, &xu);

		// update particle position (with wrap)
		xnew[id_i] = cwrap(x[id_i] + v*xu, L);
		ynew[id_i] = cwrap(y[id_i] + v*yu, L);
	}
	if(connectedSets != 0)
		connectedSets->squash();
}

void vicsekMetricThreaded(
	// model parameters
	const   size_t N,            // number of particles
	const   double L,            // linear size
	const   double v,            // particle velocity
	const   double eta,          // noise parameter (stdandard deviation)
	// variables/buffers
	FlockState & fs,
	double* const		   hAgg, // Average of neighbours for each particle
	// random number generation
	std::mt19937_64& rng,            // 64-bit Mersenne Twister
	const   int updateMethod,  //0 = backwards, 1 = forwards
	VpTreeEuclidean2D & tree,
	std::vector<int>* interactionSet,
	const int numThreads
)
{
	std::uniform_real_distribution<double> dis;
	auto h = fs.h.reader(), x = fs.x.reader(), y = fs.y.reader();
	auto hnew = fs.h.writer(), xnew = fs.x.writer(), ynew = fs.y.writer();
	std::vector<VPoint<2> > points(N);
	// precalculate unit vector velocities for every particle
	std::vector<double> xangles(N);
	std::vector<double> yangles(N);
	for(size_t i = 0; i < N; ++i) {
		sincos(h[i], &yangles[i], &xangles[i]);
		points[i].values[0] = x[i];
		points[i].values[1] = y[i];
		points[i].idx = (int)i;
	}
	tree.create(points);

	struct ThreadContext
	{
		std::vector<VPoint<2> > results;
		std::vector<int> workingInteractionSet;
#ifdef DECIMATE_INLINE
		std::mt19937_64 rng;
		std::uniform_real_distribution<double> dis;
#endif
	};

	std::vector<ThreadContext> threadContexts(numThreads);
#ifdef DECIMATE_INLINE
	// Copy the state so we don't pollute the original, so that the rest of the program
	// is deterministics regardless of whether decimation is used
	auto rng2 = rng; 
#endif

	for(auto &tc : threadContexts) {
		tc.results.reserve(50);
		if(interactionSet != 0)
			tc.workingInteractionSet.reserve(100000);
#ifdef DECIMATE_INLINE
		tc.rng.seed(rng2());
#endif
	}
#pragma omp parallel for
	for(size_t i = 0; i < N; ++i) { // for each particle
		int const threadNum = omp_get_thread_num();
		auto &context = threadContexts[threadNum];
		int id_i = points[i].idx;
		context.results.clear();
		//Work out neighbours and dists. Dnc = Distance to NCth neighbour
		int num = tree.fr_search(points[i], 1.0, &context.results);

		double xanv = 0.0, yanv = 0.0;
		for(int j = 0; j < num; ++j) {
			int id_j = context.results[j].idx;
			xanv += xangles[id_j];
			yanv += yangles[id_j];

			//add to interaction set
			if(interactionSet != 0 && id_i != id_j
#ifdef DECIMATE_INLINE
			   && context.dis(context.rng) <= 0.1
#endif
			   ) {
				context.workingInteractionSet.push_back(id_i);
				context.workingInteractionSet.push_back(id_j);
			}
		}
		// Store the aggregate of all interacting particles, without noise, for use with GTE
		//hAgg[id_i] = awrap(atan2(yanv - yangles[id_i], xanv - xangles[id_i]));
		//I think these referring to id_i instead of i is causing issues with parallelisation
		//or it's the rng
		hAgg[id_i] = awrap(atan2(yanv, xanv));

		hnew[id_i] = awrap(jitter(xanv, yanv, eta, rng, dis));

		double xu, yu;
		if(updateMethod == 1) //Forwards
			sincos(hnew[id_i], &yu, &xu); // u is unit vector in direction hnew[i]
		else //Backwards
			sincos(h[id_i], &yu, &xu);

		// update particle position (with wrap)
		xnew[id_i] = cwrap(x[id_i] + v * xu, L);
		ynew[id_i] = cwrap(y[id_i] + v * yu, L);
	}
		
	if(interactionSet != 0) {
		size_t totalCount = 0;
		for(auto &tc : threadContexts) {
			totalCount += tc.workingInteractionSet.size();
		}
		interactionSet->reserve(totalCount);
		for(auto &tc : threadContexts) {
			interactionSet->insert(interactionSet->end(), tc.workingInteractionSet.begin(), tc.workingInteractionSet.end());
		}
	}
}

void vicsekTopological(
	// model parameters
	const   size_t N,            // number of particles
	const   double L,            // linear size
	const   double v,            // particle velocity
	const   double eta,          // noise parameter (stdandard deviation)
	const	double dt,			 // timestep
	// variables/buffers
	FlockState& fs,
	double* const		   hAgg, // Average of neighbours for each particle
	// random number generation
	std::mt19937_64& rng,            // 64-bit Mersenne Twister
	const	size_t NC,			 //Number of topological neighbours
	VpTreeEuclidean2D & tree,
	const   int updateMethod,  //0 = backwards, 1 = forwards
	std::vector<int>* interactionSet,
	unsigned long int* interactionCount
	)
{
	std::uniform_real_distribution<double> dis;
	auto h = fs.h.reader(), x = fs.x.reader(), y = fs.y.reader();
	auto hnew = fs.h.writer(), xnew = fs.x.writer(), ynew = fs.y.writer();
	const int NC1 = (int)NC + 1;
	std::vector<VPoint<2> > points(N);
	// precalculate unit vector velocities for every particle
	std::vector<double> xangles(N);
	std::vector<double> yangles(N);
	for(size_t i = 0; i < N; ++i) {
		sincos(h[i], &yangles[i], &xangles[i]);
		points[i].values[0] = x[i];
		points[i].values[1] = y[i];
		points[i].idx = (int)i;
	}

	tree.create(points);

	std::vector<VPoint<2> > results;
	std::vector<int> indices;
	std::vector<double> distances;
	
	for(size_t i = 0; i<N; ++i) { // for each particle
		int id_i = points[i].idx;
		results.clear();
		indices.clear();
		distances.clear();
		//Work out neighbours and dists. Dnc = Distance to NCth neighbour
		tree.knn_search(points[i], NC1, &results, &indices, &distances);

		double xanv = 0.0, yanv = 0.0;
		for(int j = 0; j < NC1; ++j) {
			int id_j = results[j].idx;
			xanv += xangles[id_j];
			yanv += yangles[id_j];

			//add to interaction set
			if(interactionSet != 0 && id_i != id_j) {
				interactionSet->push_back(id_i);
				interactionSet->push_back(id_j);
			}
			if(interactionCount != 0 && id_i != id_j)
				(*interactionCount) += 2;
		}

		// Store the aggregate of all interacting particles, without noise, for use with GTE
		hAgg[id_i] = awrap(atan2(yanv, xanv));

		hnew[id_i] = jitter(xanv, yanv, eta, rng, dis);

		double xu, yu;
		if(updateMethod == 1) //Forwards
			sincos(hnew[id_i], &yu, &xu); // u is unit vector in direction hnew[i]
		else //Backwards
			sincos(h[id_i], &yu, &xu);

		// update particle position (with wrap)
		xnew[id_i] = cwrap(x[id_i] + v*dt*xu, L);
		ynew[id_i] = cwrap(y[id_i] + v*dt*yu, L);
	}
}

void xyModel(const size_t N, const double L, const double eta, FlockState& fs, double* const hAgg, std::mt19937_64& rng)
{
	std::uniform_real_distribution<double> dis;
	auto h = fs.h.reader();
	auto hnew = fs.h.writer();

	std::vector<double> xangles(N);
	std::vector<double> yangles(N);
	for (size_t i = 0; i < N; ++i) {
		sincos(h[i], &yangles[i], &xangles[i]);
	}

	const int x_offsets[] = { 0, 1, 0, -1 };
	const int y_offsets[] = { -1, 0, 1, 0 };
	size_t L2 = static_cast<size_t>(L);
	auto calc_new_coord = [L2](size_t coord, int offset) {
		size_t ret = coord;
		if (offset < 0) {
			ret = wrap_minus(coord, L2);
		}
		else if (offset > 0) {
			ret = wrap_plus(coord, L2);
		}
		return ret;
	};
	for (size_t i = 0; i < N; ++i) {
		size_t x_coord = i % L2;
		size_t y_coord = i / L2;

		double xanv = 0.0, yanv = 0.0;
		for (size_t j = 0; j < 4; ++j) {
			size_t x1 = calc_new_coord(x_coord, x_offsets[j]);
			size_t y1 = calc_new_coord(y_coord, y_offsets[j]);
			size_t n_idx = y1 * L2 + x1;
			xanv += xangles[n_idx];
			yanv += yangles[n_idx];
		}

		hAgg[i] = awrap(atan2(yanv, xanv));
		hnew[i] = jitter(xanv, yanv, eta, rng, dis);
	}
}

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
	const int updateMethod //0 - Metric, 1, forward
	)
{
	std::uniform_real_distribution<double> dis;
	// calculate new angles and positions

	// accumulated neighbour vectors for each particle
	/*double xanv1[N];
	double yanv1[N];
	memset(xanv1, 0, N*sizeof(double));
	memset(yanv1, 0, N*sizeof(double));*/

	// precalculate unit vector velocities for every particle
	std::vector<double> xangles(N);
	std::vector<double> yangles(N);
	for(size_t i = 0; i < N; ++i) {
		sincos(h[i], &yangles[i], &xangles[i]);
	}

	for (size_t i=0;i<N;++i) { // for each particle
		double xanv = 0.0;
		double yanv = 0.0;
		//int neighbours = 0;
		// calculate average neighbour velocity for particle i
		for (size_t j=0;j<N;++j) { // for each other particle (including self)
			if(dsqwrap(x[i] - x[j], y[i] - y[j], L) < 1.0) { // if j a neighbour of i
				//xanv1[i] += xangles[j];
				xanv += xangles[j];
				//yanv1[i] += yangles[j];
				yanv += yangles[j];
				// Save result for when updating particle j
				/*if(i != j)
				{
					xanv1[j] += xangles[i];
					yanv1[j] += yangles[i];
				}*/
			}
		}
		//fmt::print("i: {} (b)\n", neighbours);

		// update angle with noise of intensity eta
		//hnew[i] = jitter(xanv1[i],yanv1[i],eta,rng); // vna is average neighbour velocity
	/*for(size_t i = 0; i < N; ++i) {
		double xanv = 0.0;
		double yanv = 0.0;
		for(size_t j = 0; j < N; ++j) {
			if(dsqwrap(x[i] - x[j], y[i] - y[j], L) < 1.0) {
				double xnv, ynv;
				sincos(h[j], &ynv, &xnv);
				xanv += xnv;
				yanv += ynv;
			}
		}*/
		hnew[i] = jitter(xanv, yanv, eta, rng, dis);
		//hnew[i] = jitter(xanv1, yanv1, eta, rng);
		//hnew[i] = awrap(atan2(yanv, xanv));

		double xu,yu;
		if(updateMethod == 1) //Forwards
			sincos(hnew[i], &yu, &xu); // u is unit vector in direction hnew[i]
		else //Backwards
			sincos(h[i], &yu, &xu);

		// update particle position (with wrap)

		xnew[i] = cwrap(x[i] + v*xu, L);
		ynew[i] = cwrap(y[i] + v*yu, L);
	}
}

void inertialSpinModel(const size_t N,
						const double L,
						const double v,
						const double sigma,
						const double chi,		// Generalised moment of inertia (=1.25)
						const double vis,		// viscosity constant (\eta in Cavagna15, renamed. (=0.3)
						const double J,			// alignment strength (=0.8)
						const double dt,
						const   size_t NC,			 // num neighbours
						FlockState& fs,
						const double * const s,
						double * const snew,
						double * const hAgg,
						std::mt19937_64& engine,	//Generator for gaussian rng
						std::normal_distribution<>& gaussian_rng, //Gaussian rng for ISM noise
						const int updateMethod,
					   VpTreeEuclidean2D& tree,
					   double & hamiltonianEnergy,
					   std::vector<int>* interactionSet,
						unsigned long int * interactionCount,
						ConnectedSets * connectedSets)
{
	auto h = fs.h.reader(), x = fs.x.reader(), y = fs.y.reader();
	auto hnew = fs.h.writer(), xnew = fs.x.writer(), ynew = fs.y.writer();
	const int NC1 = (int)NC + 1;
	std::vector<VPoint<2> > points(N);
	std::vector<double> xangles(N);
	std::vector<double> yangles(N);
	for(size_t i = 0; i < N; ++i) {
		sincos(h[i], &yangles[i], &xangles[i]);
		points[i].values[0] = x[i];
		points[i].values[1] = y[i];
		points[i].idx = (int)i;
	}

	tree.create(points);

	std::vector<VPoint<2> > results;
	std::vector<int> indices;
	std::vector<double> distances;

	//std::random_device r;
	//std::mt19937_64 engine(r());
	//std::normal_distribution<> normal_rng(0.0, 1.0);
	auto randn = [&gaussian_rng, &engine]() { return gaussian_rng(engine); };

	hamiltonianEnergy = 0;
	double inv2Chi = 1. / (2 * chi);
	for(size_t i = 0; i < N; ++i) { // for each particle
		int id_i = points[i].idx;
		results.clear();
		indices.clear();
		distances.clear();
		//Work out neighbours and dists. Dnc = Distance to NCth neighbour
		tree.knn_search(points[i], NC1, &results, &indices, &distances);

		double xanv = 0.0, yanv = 0.0;
		double dot = 0.;
		for(int j = 0; j < NC1; ++j) {
			int id = results[j].idx;
			xanv += xangles[id];
			yanv += yangles[id];

			dot += (v*xangles[id] *v* xangles[id_i] + v*yangles[id] *v* yangles[id_i]);
			//add to interaction set
			if(interactionSet != 0 && id_i != id) {
				interactionSet->push_back(id_i);
				interactionSet->push_back(id);
			}
			if(interactionCount != 0 && id_i != id)
				(*interactionCount) += 2;

			if(connectedSets != 0 && id_i != id) {
				connectedSets->connect(id_i, id);
			}
		}
		if(connectedSets != 0 && NC1 <= 1)
			connectedSets->add(id_i);

		dot *= -J / (2 * v*v);
		hamiltonianEnergy += dot + s[id_i] * s[id_i] * inv2Chi;

		// Store the aggregate of all interacting particles, without noise, for use with GTE
		hAgg[id_i] = awrap(atan2(yanv, xanv));


		double noise_x = randn() * sigma;
		double noise_y = randn() * sigma;


		double sxv_x, sxv_y;
		crossProduct_sxv(s[id_i], xangles[id_i] * v, yangles[id_i] * v, sxv_x, sxv_y);
		double new_x = xangles[id_i] * v + (dt / chi) * sxv_x;
		double new_y = yangles[id_i] * v + (dt / chi) * sxv_y;

		//Velocity update
		hnew[id_i] = awrap(atan2(new_y, new_x));


		//or do it via matrix rotations:
		/*double spin_angle = s[id_i] * dt / chi;
		double tmp_x = xangles[id_i] * v;
		double tmp_y = yangles[id_i] * v;
		double new_x = tmp_x * cos(spin_angle) - tmp_y * sin(spin_angle);
		double new_y = tmp_x * sin(spin_angle) + tmp_y * cos(spin_angle);
		hnew[id_i] = awrap(atan2(new_y, new_x));*/

		//Cross product of v_i with consensus vector
		double hxa = crossProduct_vxv(xangles[id_i] * v, yangles[id_i] * v, xanv*v, yanv*v);
		//Cross product of v_i with noise vector
		double hxn = crossProduct_vxv(xangles[id_i] * v, yangles[id_i] * v, noise_x, noise_y);
		//Spin update
		snew[id_i] = (1.0 - vis * dt / chi) * s[id_i] + (J*dt / (v*v)) * hxa + (std::sqrt(dt) / v)*hxn;
		
		double xu, yu;
		if(updateMethod == 1) //Forwards
			sincos(hnew[id_i], &yu, &xu); // u is unit vector in direction hnew[i]
		else //Backwards
			sincos(h[id_i], &yu, &xu);
		
		// update particle position (with wrap)
		//Disable to make "on-lattice"
		xnew[id_i] = cwrap(x[id_i] + v*dt*xu, L);
		ynew[id_i] = cwrap(y[id_i] + v*dt*yu, L);
		//xnew[id_i] = x[id_i];
		//ynew[id_i] = y[id_i];
	}
	if(connectedSets != 0)
		connectedSets->squash();
}

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
	)
{
	std::uniform_real_distribution<double> dis;
	// calculate new angles and positions

	for(size_t i = 0; i<N; ++i) { // for each particle

		// calculate average neighbour velocity

		double xanv = 0.0;
		double yanv = 0.0;
		for(size_t j = 0; j<N; ++j) { // for each other particle (including self)
			if(dsqwrap(x[i] - x[j], y[i] - y[j], L) < 1.0) { // if j a neighbour of i
				double xnv, ynv;
				sincos(h[j], &ynv, &xnv); // nv is unit vector in direction h[j]
				xanv += xnv;
				yanv += ynv;
			}
		}

		// update angle with noise of intensity eta

		hnew[i] = jitter(xanv, yanv, eta, rng, dis); // vna is average neighbour velocity

//		double xu, yu;
//#ifdef VICSEK_UPDATE_FORWARDS
//		sincos(hnew[i], &yu, &xu); // u is unit vector in direction hnew[i]
//#else 
//#ifdef VICSEK_UPDATE_BACKWARDS
//		sincos(h[i], &yu, &xu);
//#endif
//#endif

		// update particle position (with wrap)

		xnew[i] = x[i];
		ynew[i] = y[i];
		//xnew[i] = cwrap(x[i] + v*xu, L);
		//ynew[i] = cwrap(y[i] + v*yu, L);
	}
}

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
	std::mt19937_64& rng            // 64-bit Mersenne Twister
)
{
	std::uniform_real_distribution<double> dis;
	// calculate new angles and positions
	for (size_t i=0;i<N;++i) { // for each particle

		// calculate average neighbour velocity

		double avel = 0.0;
		for (size_t j=0;j<N;++j) { // for each other particle (including self)
			if (dsqwrap(x[i]-x[j],y[i]-y[j],L) < 1.0) { // if j a neighbour of i
				avel += sin(h[j]-h[i]); // angular velocity term
			}
		}

		// update angle with noise of intensity eta

		hnew[i] = awrap(h[i] + w*avel + eta*(dis(rng)-0.5));

		double xu,yu;
		sincos(hnew[i],&yu,&xu); // u is unit vector in direction hnew[i]

		// update particle position (with wrap)

		xnew[i] = cwrap(x[i] + v*xu,L);
		ynew[i] = cwrap(y[i] + v*yu,L);
	}
}

void rotateModel(const size_t N, const double L, const double rotation, const double * const h, const double * const x, const double * const y,
				 double * const hnew, double * const xnew, double * const ynew)
{
	double sinTheta, cosTheta;
	//Rotation factors for 2D rot matrix
	sincos(rotation, &sinTheta, &cosTheta);
	//Translation component of 3D homogeneous rot matrix, to rotate about center of L-box, rather than origin
	double xTrans = (-L * 0.5)*cosTheta + (L*0.5)*sinTheta + L * 0.5;
	double yTrans = (-L * 0.5)*sinTheta - (L*0.5)*cosTheta + L * 0.5;

	for(size_t i = 0; i < N; ++i) {
		xnew[i] = cwrap(x[i] * cosTheta - y[i] * sinTheta + xTrans, L);
		ynew[i] = cwrap(x[i] * sinTheta + y[i] * cosTheta + yTrans, L);
		hnew[i] = awrap(h[i] + rotation);
	}
}

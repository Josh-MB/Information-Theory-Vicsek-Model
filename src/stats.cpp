#include "../include/stats.hpp"
#include "../include/nearestNeighbour.hpp"
#include "../include/vpUtil.hpp"
#include "../include/utils.hpp"

#include <vector>
#include <stdlib.h>

#ifdef ___OSX___
	#define sincos __sincos
#endif

void order_param(
	const   size_t        N,      // number of particles
	const   double* const h,      // angles, in range [-pi,pi]
	double* const         zx,     // order parmeter x coord
	double* const         zy,      // order parmeter y coord
	double  const		  cumulRot //Cumulative frame rotation
)
{
	// calculate order parameter = mean particle velocity vector
	// use atan2(zy,zx) for angle in [-pi,pi]
	// use hypot(zx,zy) for magnitude

	*zx = 0.0;
	*zy = 0.0;
	for (size_t i=0;i<N;++i) { // for each particle
		double ux,uy;
		sincos(h[i]-cumulRot,&uy,&ux); // u is unit vector in direction h[i]
		*zx += ux;
		*zy += uy;
	}
	*zx /= (double)N;
	*zy /= (double)N;
}

void nn_histogram(
	const   size_t        N,     // number of particles
	const   double        L,     // linear size
	const   double* const x,     // x coords, in range [0,L]
	const   double* const y,     // y coords, in range [0,L]
	uint64_t*  const         count  // number of neighbours histogram (accumulate)
)
{
	for (size_t i=0;i<N;++i) { // for each particle
		uint64_t n = 0;
		for (size_t j=0;j<N;++j) { // for each other particle in pair
			if (j == i) continue;  // excluding self
			//Boolean expressions equate to 0 or 1, so just add
			//result to n, rather than only adding 1 if true
			n+=(dsqwrap(x[i]-x[j],y[i]-y[j],L) < 1.0);
		}
		++count[n];
	}
}

void MI_histogram(
	const  size_t        N,      // number of particles
	const  double* const h,      // angles, in range [-pi,pi]
	const  double        cumulativeRotation,
	const  size_t        B,      // number of bins
	uint64_t* const         count,  // 2-dim histogram (accumulate)
	size_t* const	     bins,
	std::vector<int> const & interactionSet
)
{
	//Precalculate all the bins, rather than as needed
	for(size_t i=0;i<N;++i)
		bins[i] = binit(awrap(h[i]-cumulativeRotation),B);
	for(size_t i = 0; i < interactionSet.size(); i+=2) {
		const size_t ib = bins[interactionSet[i]];
		const size_t jb = bins[interactionSet[i + 1]];
		++count[ib + B*jb];
	}
	/*for (size_t i=0;i<N;++i) { // for each particle
		//const size_t ib = bins[i];//binit(h[i],B);
		for (size_t j=0;j<N;++j) { // for each other particle in pair
			if (j == i) continue;  // excluding self
			//if (dsqwrap(x[i]-x[j],y[i]-y[j],L) < 1.0) { // if j a neighbour of i
			//	const size_t jb = bins[j];//binit(h[j],B);
				++count[binit2(h[i]-h[j], B)];
			//}
		}
	}*/
}

void MI_1D_histogram(
	const  size_t        N,      // number of particles
	const  double        L,      // linear size
	const  double* const h,      // angles, in range [-pi,pi]
	const  double* const x,      // x coords, in range [0,L]
	const  double* const y,      // y coords, in range [0,L]
	const  size_t        B,      // number of bins
	uint64_t* const         count,  // 1-dim histogram (accumulate)
	uint64_t* const		 marginalCount
)
{
	for (size_t i=0;i<N;++i) { // for each particle
		++marginalCount[binit(h[i], B)];
		for (size_t j=0;j<N;++j) { // for each other particle in pair
			if (j == i) continue;  // excluding self
			if(dsqwrap(x[i] - x[j], y[i] - y[j], L) < 1.0) {
				double diff = std::abs(h[i] - h[j]);
				++count[binit3(std::min(diff, TWOPI - diff), B)];
			}
			
		}
	}
}

void GTE_1D_histogram(
	const  size_t        N,      // number of particles
	const  double* const h,      // angles, in range [-pi,pi]
	const  double* const hnew,
	const  size_t        B,      // number of bins
	uint64_t* const         count,  // 1-dim histogram (accumulate)
	const  double		 rotation
)
{
	for(size_t i = 0; i < N; ++i) {
		double diff = std::abs(awrap(hnew[i] - rotation - h[i]));
		++count[binit3(std::min(diff, TWOPI - diff), B)];
	}
}

void GTE_2D_histogram(
	const  size_t        N,      // number of particles
	const  double* const h,      // angles, in range [-pi,pi]
	const  double* const hnew,
	const  size_t        B,      // number of bins
	uint64_t* const         count,  // 2-dim histogram (accumulate)
	const  double		 rotation
)
{
	for(size_t i = 0; i < N; ++i) {
		++count[binit(hnew[i] - rotation, B)*B+binit(h[i], B)];
	}
}

void MI_histogram_shuffled(
	const  size_t        N,      // number of particles
	const  size_t 	     U,      // number of updates
	const  double* const h,      // angles, in range [-pi,pi]
	const  size_t        B,      // number of bins
	uint64_t* const         count,  // 2-dim histogram (accumulate)
	size_t* const	     bins,
	std::vector<std::vector<int> > const & fullInteractionSet,
	std::mt19937_64& rng
)
{
	std::vector<size_t> jb;
	std::vector<size_t> ib;

	for(size_t u = 0; u < U; ++u) {
		//Precalculate all the bins, rather than as needed
		for(size_t i=0;i<N;++i)
			bins[i] = binit(h[u+i*U],B);
	
		for(size_t i = 0; i < fullInteractionSet[u].size(); i += 2) {
			ib.push_back(bins[fullInteractionSet[u][i]]);
			jb.push_back(bins[fullInteractionSet[u][i + 1]]);
		}
	}

	//Shuffle y
	std::shuffle(jb.begin(), jb.end(), rng);

	//Construct histogram
	for(size_t i = 0; i < ib.size(); ++i) {
		++count[ib[i]+B*jb[i]];
	}
}

void MI_histogram_all_pairs(
	const  size_t        N,      // number of particles
	const  double* const h,      // angles, in range [-pi,pi]
	const  size_t        B,      // number of bins
	uint64_t* const         count,  // 2-dim histogram (accumulate)
	size_t* const	     bins
)
{
	//Precalculate all the bins, rather than as needed
	for(size_t i = 0; i<N; ++i)
		bins[i] = binit(h[i], B);
	for(size_t i = 0; i<N; ++i) { // for each particle
		const size_t ib = bins[i];//binit(h[i],B);
		for(size_t j = 0; j<N; ++j) { // for each other particle in pair
			if(j == i) continue;  // excluding self
								  //if (dsqwrap(x[i]-x[j],y[i]-y[j],L) < 1.0) { // if j a neighbour of i
			const size_t jb = bins[j];//binit(h[j],B);
			++count[ib + B*jb];
			//}
		}
	}
}

void TE_histogram(
	const  size_t        N,      // number of particles
	const  double* const h,      // angles, in range [-pi,pi]
	const  double* const hnew,   // new angles
	const  double        rotation,
	const  double        cumulativeRotation,
	const  size_t        B,      // number of bins
	uint64_t* const         count,  // 3-dim histogram (accumulate)
	size_t* const	     bins,
	std::vector<int> const & interactionSet
)
{
	//Precalculate all the bins, rather than as needed
	for(size_t i=0;i<N;++i)
		bins[i] = binit(awrap(h[i]-cumulativeRotation),B);

	for(size_t i = 0; i < interactionSet.size(); i += 2) {
		const size_t ib = bins[interactionSet[i]];
		const size_t jb = bins[interactionSet[i + 1]];
		const size_t inewb = binit(awrap(hnew[interactionSet[i]] - rotation - cumulativeRotation), B);
		++count[inewb + B*(ib + B*jb)];
	}
}

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
)
{
	std::vector<size_t> ib[3];
	
	for(size_t u=0; u < U-1; ++u) {
		for(size_t i=0;i<N;++i)
			bins[i] = binit(h[u+i*U],B);

		for(size_t i = 0; i < fullInteractionSet[u].size(); i += 2) {
			ib[0].push_back(bins[fullInteractionSet[u][i]]);
			ib[1].push_back(bins[fullInteractionSet[u][i + 1]]);
			ib[2].push_back(binit(h[(u + 1)+(fullInteractionSet[u][i]) * U], B));
		}
	}

	std::shuffle(ib[d].begin(), ib[d].end(), rng);

	for(size_t i=0; i < ib[0].size(); ++i) {
		++count[ib[2][i]+B*(ib[0][i]+B*ib[1][i])];
	}
}

void insort(
	const  size_t   N,            // number of elements to sort
	double*         const elem,   // elements to sort
	size_t*         const idx     // indices of elements
)
{
	// simple insertion sort

	for (size_t j=0;j<N;++j) idx[j] = j;
	for (size_t j=0;j<N;++j) {
		const double d = elem[j];
		size_t k = j;
		while (k>0 && elem[k-1]>d) {
			elem[k] = elem[k-1];
			idx[k]  = idx[k-1];
			--k;
		}
		elem[k] = d;
		idx[k]  = j;
	}
}

void siftDown(
	double*		const elem,
	size_t*		const idx,
	int 		root,
	int		bottom
)
{
	//Helper function to heapsort algorithm
	int maxChild = root * 2 + 1;

	if(maxChild < bottom)
	{
		int otherChild = maxChild + 1;
		maxChild = (elem[otherChild] > elem[maxChild])?otherChild:maxChild;
	}
	else
	{
		if(maxChild > bottom)
			return;
	}
	
	if(elem[root] >= elem[maxChild])
		return;

	double tempElem = elem[root];
	elem[root] = elem[maxChild];
	elem[maxChild] = tempElem;

	size_t tempIdx = idx[root];
	idx[root] = idx[maxChild];
	idx[maxChild] = tempIdx;

	siftDown(elem, idx, maxChild, bottom);
}

void heapsort(
	const size_t	N,
	double*		const elem,
	size_t*		const idx
)
{
	int i;
	size_t tempIdx;
	double tempElem;	
	
	for(i = (int)N / 2; i >= 0; --i)
		siftDown(elem, idx, i, (int)N - 1);

	for(i = (int)N - 1; i >= 1; --i)
	{
		//Swap
		tempElem = elem[0];
		elem[0] = elem[i];
		elem[i] = tempElem;
		tempIdx = idx[0];
		idx[0] = idx[i];
		idx[i] = tempIdx;
	
		siftDown(elem, idx, 0, i-1);
	}
}

void GTE_histogram(
	const   size_t        N,      // number of particles
	const   double        L,      // linear size
	const   double* const h,      // angles, in range [-pi,pi]
	const   double* const hnew,   // new angles
	const   double* const x,      // x coords, in range [0,L]
	const   double* const y,      // y coords, in range [0,L]
	const   size_t        K,      // number of nearest neighbours in source (including self!)
	double*         const dsq,    // square distances buffer
	size_t*         const idx,    // indices buffer
	const   size_t        B,      // number of bins
	uint64_t*  const         count, // (K+1)-dim histogram (accumulate)
	size_t* const	     bins
)
{
	//Precalculate all the bins, rather than as needed
	for(size_t i=0;i<N;++i)
		bins[i] = binit(h[i],B);
	for (size_t i=0;i<N;++i) { // for each particle

		// simple insertion sort on square distances (maybe qsort better here...)
		// note: i itself will sort to first in array,
		// that is, idx[0] = i,  since distance = 0 !

		for (size_t j=0;j<N;++j) dsq[j] = dsqwrap(x[i]-x[j],y[i]-y[j],L);
		insort(N,dsq,idx);
		//Outputs differ. Maybe due to instability of heapsort
		//as opposed to stability of insort.
		//Doesn't matter, won't be sorting later anyway
		//heapsort(N,dsq,idx);

		// accumulate count in appropriate bin
		// note 1: dimension of histogram = B^(K+1)
		// note 2: first array dimension is hnew[i], second is h[i], subsequent are K-1 nearest neighbours

		size_t Bk = 1;
		size_t ib = binit(hnew[i],B);
		for (size_t k=0;k<K;++k) { // for K nearest neighbours (including self)
			Bk *= B;
			const size_t jk = bins[idx[k]];//binit(h[idx[k]],B); // k-th nearest neighbour
			ib += Bk*jk;
		}
		++count[ib];
	}
}

void GTE_avg_histogram(
	const  size_t        N,      // number of particles
	const  double* const h,      // angles, in range [-pi,pi]
	const  double* const hnew,   // new angles
	const  double	     rotation,
	const  double        cumulativeRotation,
	const  double* const hAgg,   // neighbour angles
	const  size_t        B,      // number of bins
	uint64_t* const         count  // 3-dim histogram (accumulate)
)
{
	//Precalculate all the bins, rather than as needed
	for(size_t i = 0; i<N; ++i) { // for each particle
		const size_t inewb = binit(awrap(hnew[i] - rotation - cumulativeRotation), B);
		const size_t ib = binit(awrap(h[i]-cumulativeRotation), B);
		const size_t jb = binit(awrap(hAgg[i]-cumulativeRotation),B);
		++count[inewb + B*(ib + B*jb)];
	}
}

void GTE_avg_histogram_shuffled(
	const  size_t        N,      // number of particles
	const  size_t        U,      // number of updates
	const  double* const h,      // angles, in range [-pi,pi]
	const  double* const hAgg,   // neighbour angles
	const  size_t        B,      // number of bins
	uint64_t* const         count,  // 3-dim histogram (accumulate)
	const  int	     d,	     // dimension to shuffle (0 = x, 1 = y, 2 = w),
	std::mt19937_64& rng	
)
{
	std::vector<uint64_t> ib[3];
	ib[0].reserve(N*U);
	ib[1].reserve(N*U);
	ib[2].reserve(N*U);
	for(size_t u=0; u < U-1; ++u) {
		for(size_t i = 0; i<N; ++i) {
			ib[0].push_back(binit(h[u+i*U], B));
			ib[1].push_back(binit(hAgg[u+i*U], B));
			ib[2].push_back(binit(h[u+1+i*U], B));
		}
	}
	std::shuffle(ib[d].begin(), ib[d].end(), rng);

	//Precalculate all the bins, rather than as needed
	for(size_t i = 0; i<ib[0].size(); ++i) { // for each particle
		++count[ib[2][i] + B*(ib[0][i] + B*ib[1][i])];
	}
}

void calcCenterMass(const size_t N, double const * const x, double const * const y, const double L,
					double &cx, double &cy)
{
	double xi[2] = { 0., 0. };
	double zeta[2] = { 0., 0. };

	for(size_t i = 0; i < N; ++i) {
		double theta = TWOPI * x[i] / L;
		xi[0] += cos(theta);
		zeta[0] += sin(theta);

		theta = TWOPI * y[i] / L;
		xi[1] += cos(theta);
		zeta[1] += sin(theta);
	}

	xi[0] /= (double)N;
	xi[1] /= (double)N;
	zeta[0] /= (double)N;
	zeta[1] /= (double)N;

	cx = L * (atan2(-zeta[0], -xi[0]) + PI) / TWOPI;
	cy = L * (atan2(-zeta[1], -xi[1]) + PI) / TWOPI;
}

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
)
{
	std::vector<double> adj_x0(N), adj_y0(N);
	std::vector<double> adj_xt(N), adj_yt(N);

	//Adjust both timesteps such that their origins are the CM
	//for those timesteps
	for(size_t i = 0; i < N; ++i) {
		adj_x0[i] = offsetwrap(cm_x0, x0[i], L);
		adj_y0[i] = offsetwrap(cm_y0, y0[i], L);
		adj_xt[i] = offsetwrap(cm_xt, xt[i], L);
		adj_yt[i] = offsetwrap(cm_yt, yt[i], L);
	}
	
	std::vector<double> vx0(N), vy0(N);
	std::vector<double> vxt(N), vyt(N);

	double flock_vx0 = 0., flock_vy0 = 0.;
	double flock_vxt = 0., flock_vyt = 0.;

	//Calculate velocity vectors and flock averages
	for(size_t i = 0; i < N; ++i) {
		sincos(h0[i], &vy0[i], &vx0[i]);
		sincos(ht[i], &vyt[i], &vxt[i]);
		vx0[i] *= v;
		vy0[i] *= v;
		vxt[i] *= v;
		vyt[i] *= v;
		flock_vx0 += vx0[i];
		flock_vy0 += vy0[i];
		flock_vxt += vxt[i];
		flock_vyt += vyt[i];
	}
	flock_vx0 /= (double)N;
	flock_vy0 /= (double)N;
	flock_vxt /= (double)N;
	flock_vyt /= (double)N;

	std::vector<double> dvx0(N), dvy0(N);
	std::vector<double> dvxt(N), dvyt(N);

	for(size_t i = 0; i < N; ++i) {
		dvx0[i] = (vx0[i] - flock_vx0) / v;
		dvy0[i] = (vy0[i] - flock_vy0) / v;
		dvxt[i] = (vxt[i] - flock_vxt) / v;
		dvyt[i] = (vyt[i] - flock_vyt) / v;
	}

	std::vector<double> coeffs(numR);
	for(size_t r = 0; r < numR; ++r) {
		coeffs[r] = 1 / ((double)N * 4 * PI * (double)r * (double)r * rho);
	}
	for(size_t i = 0; i < N; ++i) {
		for(size_t j = 0; j < N; ++j) {
			double dist = std::sqrt(dsqwrap(adj_x0[i] - adj_xt[j], adj_y0[i] - adj_yt[j], L));
			for(size_t r = 0; r < numR; ++r) {
				if(rVals[r] < dist && dist < rVals[r] + dr) {
					Crt_accum[r] += coeffs[r] * (dvx0[i] * dvxt[j] + dvy0[i] * dvyt[j]);
				}
			}
		}
	}
}

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
)
{
	std::vector<double> adj_x0(N), adj_y0(N);
	std::vector<double> adj_xt(N), adj_yt(N);

	//Adjust both timesteps such that their origins are the CM
	//for those timesteps
	for(size_t i = 0; i < N; ++i) {
		adj_x0[i] = offsetwrap(cm_x0, x0[i], L);
		adj_y0[i] = offsetwrap(cm_y0, y0[i], L);
		adj_xt[i] = offsetwrap(cm_xt, xt[i], L);
		adj_yt[i] = offsetwrap(cm_yt, yt[i], L);
	}

	std::vector<double> vx0(N), vy0(N);
	std::vector<double> vxt(N), vyt(N);

	double flock_vx0 = 0., flock_vy0 = 0.;
	double flock_vxt = 0., flock_vyt = 0.;

	//Calculate velocity vectors and flock averages
	for(size_t i = 0; i < N; ++i) {
		sincos(h0[i], &vy0[i], &vx0[i]);
		sincos(ht[i], &vyt[i], &vxt[i]);
		vx0[i] *= v;
		vy0[i] *= v;
		vxt[i] *= v;
		vyt[i] *= v;
		flock_vx0 += vx0[i];
		flock_vy0 += vy0[i];
		flock_vxt += vxt[i];
		flock_vyt += vyt[i];
	}
	flock_vx0 /= (double)N;
	flock_vy0 /= (double)N;
	flock_vxt /= (double)N;
	flock_vyt /= (double)N;

	std::vector<double> dvx0(N), dvy0(N);
	std::vector<double> dvxt(N), dvyt(N);

	for(size_t i = 0; i < N; ++i) {
		dvx0[i] = (vx0[i] - flock_vx0) / v;
		dvy0[i] = (vy0[i] - flock_vy0) / v;
		dvxt[i] = (vxt[i] - flock_vxt) / v;
		dvyt[i] = (vyt[i] - flock_vyt) / v;
	}

	for(size_t i = 0; i < N; ++i) {
		for(size_t j = 0; j < N; ++j) {
			if(i == j) continue;
			double dist = std::sqrt(dsqwrap(adj_x0[i] - adj_xt[j], adj_y0[i] - adj_yt[j], L));
			for(size_t k = 0; k < numK; ++k) {
				Crt_accum[k] += (1.0/(double)N) * (sin(kVals[k]*dist) / kVals[k] *dist) * (dvx0[i] * dvxt[j] + dvy0[i] * dvyt[j]);
			}
		}
	}
}

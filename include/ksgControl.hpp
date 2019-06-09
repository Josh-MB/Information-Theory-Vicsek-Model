#ifndef KSGCONTROL_HPP
#define KSGCONTROL_HPP

#include "vpUtil.hpp"
#include "defs.hpp"
#include "utils.hpp"
#include "stats.hpp"

#include <vector>
#include <stdlib.h>

/**
 * Calculates information theory metrics over a single time-step
 * \param gp - Generator params used
 * \param sp - Simulation params used
 * \param interactions - Tracked interactions for this timestep, where
 *                       (2i, 2i+1) are the two particles involved in the
 *                       ith interaction.
 * \param u - Current timestep
 * \param fs - Current flockstate
 * \param hAgg - Consensus vectors for each particle
 * \param fullGTENNidx - Buffer used for the calculating nearest K neighbours. Must be
 *                       N*K in size
 * \param posPts - Buffer used for nearest neighbour calculation. Must be N in size
 * \param itStats - In/out parameter, tracking min/max/avg values for info theory metrics
 */
void fillAndRunPairwiseKSG(GeneratorParams const& gp,
						   SimParams const& sp,
						   std::vector<int> const& interactions,
						   size_t const u,
						   FlockState & fs,
						   double const*const hAgg,
						   int * fullGTENNidx,
						   std::vector<VPoint<2> > & posPts,
						   ITStats& itStats
);

/**
 * Helper functions for calculating instantaneous metrics
 */
double fillAndRunPairwiseMI(std::vector<int> const& interactions, double const* const h, size_t const KSG_neighbours);
double fillAndRunPairwiseTE(std::vector<int> const& interactions, double const* const h, double const* const hnew, size_t const KSG_neighbours);
double fillAndRunPairwiseGTE(GeneratorParams const& gp,
	SimParams const& sp,
	std::vector<VPoint<2> >& posPts,
	double const* const x,
	double const* const y,
	double const* const h,
	double const* const hnew,
	double const* const hAgg,
	int* fullGTENNidx);

/**
 * Calculates information theory metrics over whole series, using subset decimation
 * \param gp - Generator params used
 * \param sp - Simulation params used
 * \param interactions - Tracked interactions for all timesteps, where
 *                       (2i, 2i+1) of the t-th vector are the two particles involved in the
 *                       ith interaction of time t.
 * \param fh - Entire history of the flock states
 * \param hAggFull - Consensus vectors for each particle over all time steps
 * \param cumulativeRotations - Cumulative rotations over each timestep. Only used if --rotate-frame used
 * \param rng
 * \return itStats - In/out parameter, tracking min/max/avg values for info theory metrics
 */
ITStats fillAndRunSubsetKSG(GeneratorParams const& gp,
						 SimParams const& sp,
						 std::vector<std::vector<int> > const& interactions,
						 FlockHistory const& fh,
						 double const*const hAggFull,
						 double const*const cumulativeRotations,
						 std::mt19937_64& rng            // 64-bit Mersenne Twister
						 );

/**
 * Calculates information theory metrics over whole series
 * \param gp - Generator params used
 * \param sp - Simulation params used
 * \param interactions - Tracked interactions for all timesteps, where
 *                       (2i, 2i+1) of the t-th vector are the two particles involved in the
 *                       ith interaction of time t.
 * \param fh - Entire history of the flock states
 * \param hAggFull - Consensus vectors for each particle over all time steps
 * \param cumulativeRotations - Cumulative rotations over each timestep. Only used if --rotate-frame used
 * \param fullGTENNidx - Buffer used for the calculating nearest K neighbours. Must be
 *                       N*K in size
 * \param posPts - Buffer used for nearest neighbour calculation. Must be N in size
 * \param rng
 * \return itStats - In/out parameter, tracking min/max/avg values for info theory metrics
 */
ITStats fillAndRunWholeSeriesKSG(GeneratorParams const& gp,
							  SimParams const& sp,
							  std::vector<std::vector<int> > const& interactions,
							  FlockHistory const& fh,
							  double const*const hAggFull,
							  double const*const cumulativeRotations,
							  int * fullGTENNidx,
							  std::vector<VPoint<2> > & posPts,
							  std::mt19937_64& rng
);

/**
 * Calculates information theory metrics over whole series, but relative to flock average
 * \param gp - Generator params used
 * \param sp - Simulation params used
 * \param interactions - Tracked interactions for all timesteps, where
 *                       (2i, 2i+1) of the t-th vector are the two particles involved in the
 *                       ith interaction of time t.
 * \param fh - Entire history of the flock states
 * \param hAggFull - Consensus vectors for each particle over all time steps
 * \param cumulativeRotations - Cumulative rotations over each timestep. Only used if --rotate-frame used
 * \param fullGTENNidx - Buffer used for the calculating nearest K neighbours. Must be
 *                       N*K in size
 * \param posPts - Buffer used for nearest neighbour calculation. Must be N in size
 * \param rng
 * \return itStats - In/out parameter, tracking min/max/avg values for info theory metrics
 */
ITStats fillAndRunLocalInteractionsKSG(GeneratorParams const& gp,
						SimParams const& sp,
						std::vector<std::vector<int> > const& interactions,
					    FlockHistory const& fh,
						double const*const hAggFull,
						std::vector<double>* avgDirections,
						int const*const flockIds,
						double const*const cumulativeRotations,
						int * fullGTENNidx,
						std::vector<VPoint<2> > & posPts,
						std::mt19937_64& rng
);

/**
 * Calculates information theory metrics between particle and flock avg direction 
 * over whole series. Currently only MI is implemented.
 * \param gp - Generator params used
 * \param sp - Simulation params used
 * \param interactions - Tracked interactions for all timesteps, where
 *                       (2i, 2i+1) of the t-th vector are the two particles involved in the
 *                       ith interaction of time t.
 * \param fh - Entire history of the flock states
 * \param hAggFull - Consensus vectors for each particle over all time steps
 * \param cumulativeRotations - Cumulative rotations over each timestep. Only used if --rotate-frame used
 * \param fullGTENNidx - Buffer used for the calculating nearest K neighbours. Must be
 *                       N*K in size
 * \param posPts - Buffer used for nearest neighbour calculation. Must be N in size
 * \param rng
 * \return itStats - In/out parameter, tracking min/max/avg values for info theory metrics
 */
ITStats fillAndRunFlockvsParticleKSG(GeneratorParams const& gp,
								  SimParams const& sp,
								  std::vector<std::vector<int> > const& interactions,
								  FlockHistory const& fh,
								  double const*const hAggFull,
								  std::vector<double>* avgDirections,
								  int const*const flockIds,
								  double const*const cumulativeRotations,
								  int * fullGTENNidx,
								  std::vector<VPoint<2> > & posPts,
								  std::mt19937_64 &rng
);



#endif //KSGCONTROL_HPP

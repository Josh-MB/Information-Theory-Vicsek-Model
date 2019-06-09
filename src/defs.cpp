#include "../include/defs.hpp"
#include "../include/vpUtil.hpp"
#include "../include/nearestNeighbour.hpp"
#include "../include/utils.hpp"

void restrictedConsensusVector(const size_t N, const double L, const double * const x, const double * const y, const double*const h, double * const hAgg, const size_t numNeighbours)
{
	if(numNeighbours < 1) {
		PEEXIT("Needs at least one neighbour");
	}
	std::vector<VPoint<2> > posPts(N);
	for(size_t i = 0; i < N; ++i) {
		posPts[i].values[0] = x[i];
		posPts[i].values[1] = y[i];
		posPts[i].idx = (int)i;
	}
	std::vector<int> gteAvgnnIdx(N*numNeighbours);
	getNearestNeighbours(posPts, (int)N, (int)numNeighbours, L, &gteAvgnnIdx[0]);

	for(size_t i = 0; i < N; ++i) {
		double avg = 0.;
		for(size_t j = 0; j < numNeighbours; ++j) {
			avg += h[gteAvgnnIdx[i*numNeighbours + j]];
		}
		avg /= (double)numNeighbours;
		hAgg[i] = avg;
	}
}

clara::detail::Parser generatorParamsParser(GeneratorParams & gp) {
	using namespace clara::detail;
	auto cli = Opt(gp.N, "num")["-N"]["--num-particles"]("Number of particles")
		| Opt(gp.rho, "density")["--rho"]("Particle interaction intensity")
		| Opt(gp.v, "velocity")["-v"]["--velocity"]("Particle velocity")
		| Opt(gp.eta, "noise")["--eta"]("Noise")
		| Opt(gp.S, "timesteps")["-S"]["--skipsteps"]("Number of updates to skip")
		| Opt(gp.seed, "number")["--seed"]("Random seed (0 for unpredictable)")
		| Opt(gp.imethod, "method")["--imethod"]("Interaction Method, 0 = metric, 1 = topological")
		| Opt(gp.umethod, "method")["--umethod"]("Update Method, 0 = backwards, 1 = forwards")
		| Opt(gp.topo_neighbours, "num")["--topo-neighbours"]("Number of neighbours for topological updating")
		| Opt(gp.rotate_frame)["--rotate-frame"]("Rotate the reference frame between updates")
		| Opt(gp.discretise, "bins")["--discretise"]("Discretise model into number of bins")
		| Opt(gp.chi, "chi")["--ism-chi"]("Generalised moment of inertia for ISM")
		| Opt(gp.J, "J")["--ism-align-strength"]("Alignment strength for ISM")
		| Opt(gp.viscosity, "viscosity")["--ism-viscosity"]("Viscosity coefficient for ISM (\\eta in their equations)")
		| Opt(gp.dt_factor, "scale")["--ism-delta-t-factor"]("Scaling factor for Delta T for ISM")
		| Opt(gp.force_dt, "time")["--force-delta-t"]("Override Delta T calculation and use this value");
	return cli;
}

clara::detail::Parser simParamsParser(SimParams & sp) {
	using namespace clara::detail;
	auto cli = Opt(sp.sims[MIBIN])["--MI"]("Calculate MI")
		| Opt(sp.sims[TEBIN])["--TE"]("Calculate TE")
		| Opt(sp.sims[GTEBIN])["--GTE"]("Calculate GTE")
		| Opt(sp.sims[PARAMS])["--params"]("Output parameters with no information theoretic calculations")
		| Opt(sp.U, "timesteps")["-U"]["--updatesteps"]("Number of update steps")
		| Opt(sp.B, "bins")["-B"]["--bins"]("Number of bins to use in histograms")
		| Opt(sp.angleWrite)["--write-angle-to-file"]("Saves entire angle history to the output file")
		| Opt(sp.seed, "number")["--seed"]("Random seed (0 for unpredictable)")
		| Opt(sp.reseed)["--update-initial-state"]("Overwrite the initialState file with the final state here, for use with future runs")
		| Opt(sp.ksg_t)["--calc-metrics-per-timestep"]("Calculate metrics per timestep")
		| Opt(sp.ksg_w)["--calc-metrics"]("Calculate metrics (uses whole series)")
		| Opt(sp.ksg_local, "mode")["--calc-metrics-special"]("Extra calculation modes. 1 = Calculate metrics relative to flock heading. 2 = Calculate metrics with \theta_J = flock heading")
		| Opt(sp.ksg_gte_dims, "dimensions")["--GTE-KSG-dimensions"]("Number of dimensions to use for KSG calculation. 0='Full' N-dimensional GTE, 2=Tgl^2D (no-consensus-vector optimisation), 3=Tgl^3D (consensus-vector optimisation)")
		| Opt(sp.hist_gte_dims, "dimensions")["--GTE-hist-dimensions"]("Number of dimensions to use for histogram estimation of GTE, 0=disabled, 1 for calculating long-term statistics or 2 for short-term statistics")
		| Opt(sp.GKavg, "number")["--GTE-consensus-average"]("How many neighbours to use when calculating the consensus vector for metrics. 0 = use all in range")
		| Opt(sp.KSG_neighbours, "number")["--KSG-neighbours"]("How many neighbours to use for the KSG estimator")
		| Opt(sp.shuffle)["--shuffle"]("Shuffle y timeseries for testing")
		| Opt(sp.te_shuffle_dim, "dimension")["--shuffle-dimension"]("Use a different dimension for shuffling, 0=x, 1=y, 2=w")
		| Opt(sp.subset_size, "number or percentage")["--subset-size"]("Size of subset to use (0 = no subsets, otherwise min 1000, or between 0 and 1 for percentage of total points)")
		| Opt(sp.subset_runs, "runs")["--subset-runs"]("How many subsets to generate and process (Must be >0 if --subset-size > 0)")
		| Opt(sp.record_T_steps, "gap")["--record-nth-step"]("Record every nth timestep for data analysis");
	return cli;
}
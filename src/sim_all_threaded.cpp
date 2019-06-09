#include "../include/psi.hpp"
#include "../include/model.hpp"
#include "../include/stats.hpp"
#include "../include/dataSet.hpp"
#include "../include/metrics.hpp"
#include "../include/logFileHelper.hpp"
#include "../include/nearestNeighbour.hpp"
#include "../include/metricTests.hpp"
#include "../include/vpUtil.hpp"
#include "../include/ksgControl.hpp"
#include "../include/version.hpp"
#include "../include/connectedFlocks.hpp"
#include "../include/vpTree.hpp"
#include "../include/utils.hpp"
#include "../include/itHistograms.hpp"

#include <fmt/format.h>
#include <clara.hpp>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <set>
#include <algorithm>
#include <random>
#include <omp.h>

int sim_all_threaded(int argc, char* argv[])
{
	std::string initialState_nonConst, ofile_nonConst = "incorrectName";
	int threadLimit_nonConst = 1;
	SimParams sp{};
	{
		using namespace clara::detail;
		bool showHelp = false;
		auto cli = Help(showHelp)
			| Opt(initialState_nonConst, "filename")["-i"]["--in-file"]("Input file containing initial state. Leave blank for random initialisation.").required()
			| Opt(ofile_nonConst, "filename")["-o"]["--out-file"]("Base name for output data. Will prepend mibin/tebin/gtebin and append .bin or .log").required()
			| Opt(threadLimit_nonConst, "num threads")["--threads"]("Max number of threads to use")
			| simParamsParser(sp);
		auto result = cli.parse(Args(argc, argv));
		if (!result) {
			fmt::print("Error in command line: {}\n", result.errorMessage());
			return EXIT_FAILURE;
		}
		if (showHelp) {
			std::cerr << cli << std::endl;
			return EXIT_SUCCESS;
		}
	}
	const std::string initialState = initialState_nonConst, ofile = ofile_nonConst;
	const int threadLimit = threadLimit_nonConst;

	FILE* logFiles[SIM_COUNT];
	FILE* binFiles[SIM_COUNT];

	int numThreads = omp_get_max_threads();
	omp_set_num_threads(numThreads);
	if(threadLimit > 0) {
		omp_set_num_threads(threadLimit);
		numThreads = threadLimit;
		fmt::print("Limiting threads used to {} (Using {}, max {})\n", threadLimit, numThreads, omp_get_max_threads());
	}
	fmt::print("Max threads: {} ({})\n", omp_get_max_threads(), numThreads);
	fmt::print("Using {} threads\n", numThreads);

	/***
	 * Open up files for writing
	 ***/
	openFiles(sp.sims, logFiles, binFiles, ofile.c_str());

	/**
	 * Load initial state and run U updates
	 */
	if(initialState.empty()) PEEXIT("No initial state to use");
	
	fmt::print("Git Head: {}\n", GIT_HEAD);
	fmt::print("Git Log: {}\n", GIT_LOG);
	fmt::print("Git Diff: {}\n", GIT_DIFF);

	fmt::print("Initialising model\n");

	sp.seed = calcSeed(sp.seed);
	std::mt19937_64 rng(sp.seed);
	fmt::print("seed : {}\n", sp.seed);

	FILE* initFile = fopen(initialState.c_str(), "rb");
	if(initFile == NULL) PEEXIT("failed to open initialisation file");

	auto binVars = readBinaryHeader(initFile);
	GeneratorParams gp = std::get<0>(binVars);

	size_t const recordedU = sp.U / sp.record_T_steps; // Number of timesteps we record data from

	fmt::print("linear size (L)    : {}\n", gp.L);
	fmt::print("recordU      : {}\n", recordedU);

	/**
	 * Set up data buffers
	 */
	FlockHistory entire_history(gp.N, recordedU);
	FlockState fs(gp.N);

	// Heading aggregates
	double* const hAgg = (double* const)calloc(gp.N, sizeof(double));
	double* const hAgg_history = (double* const)calloc(gp.N*recordedU, sizeof(double)); // note: calloc initialises to zero - this is what we want!

	// rotations[t] is the alpha term in the equation
	// \theta_i^{t+1} = <\theta_t>_{k_nn} + noise + \alpha 
	double* const rotations = (double* const)calloc(recordedU, sizeof(double));
	double* const cumulativeRotations = (double* const)calloc(recordedU + 1, sizeof(double));

	int* fullgtennIdx = new int[gp.N*GK]; //Buffer used for NN neighbour calculation
	std::vector<VPoint<2> > posPts(gp.N); //Buffer used for NN neighbour calculation

	// Buffers used for record interacting particles over all timesteps
	std::vector<std::vector<int> > allInteractions(recordedU);
	for (auto& v : allInteractions)
		v.reserve(gp.N*gp.topo_neighbours);

	std::vector<int> oneTimeInteractions; // Dummy buffer for recording particles over a single timestep

	// Histogram accumulation buffers
	ITHistograms itHists(gp, sp);

	// Reusable VP Tree
	VpTreeEuclidean2D vicsekTree(gp.L);

	// Buffers for calculating connected sub-flocks and their average directions
	ConnectedSets connectedSet(gp.N);
	int* const flockIds = (int* const)calloc(recordedU * gp.N, sizeof(int));
	std::vector<std::vector<double>> avgDirections(recordedU);

	// Buffer for calculating info theory stats every timestep
	ITStats per_timestep_stats;

	// Buffer for calculating 1D GTE using KL estimator
	DataSet<1> KL1Ddata((int)(recordedU * gp.N), (int)sp.KSG_neighbours);
	int klIdx = 0;

	readBinaryContents(initFile, fs);
	fclose(initFile);

	fs.swapBuffers();
	duplicate_reader_to_writer(fs.h);
	duplicate_reader_to_writer(fs.x);
	duplicate_reader_to_writer(fs.y);

	/***
	* Write out header details in log and bin files
	***/
	writeFileHeaders(logFiles, binFiles, gp, sp, ofile.c_str());

	/*********************
	* Perform  analysis steps
	*********************/
	fmt::print("Performing analysis updates\n");

	/**********************************
	 *
	 * Run model
	 *
	 **********************************/
	cumulativeRotations[0] = 0;
	size_t buffer_u = 0; // Index to buffers for this u. buffer_u = u / record_T_steps
	for(size_t u = 0; u < sp.U; ++u) { // for each update
		std::vector<int>* interactionBuffer = &oneTimeInteractions;
		auto h = fs.h.reader(), x = fs.x.reader(), y = fs.y.reader();
		auto hnew = fs.h.writer();

		if(u % sp.record_T_steps == 0) {
			if(sp.sims[MIBIN] || sp.sims[TEBIN])
				interactionBuffer = &(allInteractions[buffer_u]);
		}
		interactionBuffer->clear();
		interactionBuffer->reserve(gp.N*gp.topo_neighbours);

		connectedSet.clear();
		//fmt::print("Timestep: {}\n", u);
		if(gp.imethod == 0) {
			if(gp.rotate_frame == 0)
				vicsekMetricThreaded(gp.N, gp.L, gp.v, gp.eta, fs, hAgg, rng, gp.umethod, vicsekTree, interactionBuffer, numThreads);
		}
		/*else if(p.imethod == 1) {
			vicsekTopological(p.N, p.L, p.v, p.eta, dt, h, x, y, hnew, xnew, ynew, hAgg, &rng, p.topo_neighbours, vicsekTree, p.umethod, interactionBuffer);
		}*/
		else
			PEEXIT("Invalid imethod");

		if(u % sp.record_T_steps == 0) {
			calcAverageFlockDirections(connectedSet, h, avgDirections[buffer_u]);
			memcpy(&(flockIds[buffer_u*gp.N]), &(connectedSet.particleMap[0]), gp.N * sizeof(int));
			cumulativeRotations[buffer_u + 1] = awrap(cumulativeRotations[buffer_u] + rotations[buffer_u]);

			if(sp.GKavg != 0) {
				restrictedConsensusVector(gp.N, gp.L, x, y, h, hAgg, sp.GKavg);
			}
			calc_and_write_order_param(gp.N, h, cumulativeRotations[buffer_u], sp.sims, binFiles);

			entire_history.recordState(fs, buffer_u);
			for(size_t i = 0; i < gp.N; ++i) {
				hAgg_history[buffer_u + i*recordedU] = hAgg[i];
			}

			fmt::print("Interactions (t={}): {}\n", (int)u, (int)interactionBuffer->size() / 2);

			itHists.accumulateHists(gp, fs, rotations[buffer_u], cumulativeRotations[buffer_u], allInteractions[buffer_u]);
			for(size_t i = 0; i < gp.N; ++i) {
				KL1Ddata.setData(awrap(hnew[i] - rotations[buffer_u] - h[i]), klIdx, 0);
				++klIdx;
			}

			// Run Metrics
			if(sp.ksg_t) {
				fillAndRunPairwiseKSG(gp, sp, allInteractions[buffer_u], u, fs, hAgg,
									  fullgtennIdx, posPts, per_timestep_stats);
			}

			buffer_u++;
		}

		// swap buffers
		fs.swapBuffers();

		progrep("sim ", u, sp.U);
	}
	double kl1DGTE = 0.0;
	if(sp.sims[GTEBIN] && !sp.shuffle && sp.hist_gte_dims == 1) {
		KL1Ddata.calcBallsizes();
		kl1DGTE = -psi(KL1Ddata.getK()) + psi(KL1Ddata.getN()) + KL1Ddata.averageBallsize() - log(gp.eta);
	}
	//If want to use this, need to fix miHist and gteHist. They're currently set up for 1D histogram, but
	//these methods assume they're multidimensional
	//if(p.shuffle) {
	//	if(p.sims[MIBIN]) {
	//		MI_histogram_shuffled(p.N, recordedU, count, p.B, miHist, bins, allInteractions, rng);
	//	}
	//	if(p.sims[TEBIN]) {
	//		TE_histogram_shuffled(p.N, recordedU, count, p.B, teHist, bins, p.te_shuffle_dim, allInteractions, rng);
	//	}
	//	if(p.sims[GTEBIN] && p.gte_avg == 1) {
	//		GTE_avg_histogram_shuffled(p.N, recordedU, count, hAgg_history, p.B, gteHist, p.te_shuffle_dim, rng);
	//	}
	//	//Todo: Need to shuffle GTE
	//}

	// Reseed the state file for cooling regime
	if(sp.reseed) {
		reseed(initialState.c_str(), gp, sp, fs);
	}

	// Free no longer used memory
	free(hAgg);

	ITStats subset_results;
	ITStats whole_results;
	ITStats localised_results;
	if(sp.subset_size != 0)
		subset_results = fillAndRunSubsetKSG(gp, sp, allInteractions, entire_history, hAgg_history, cumulativeRotations, rng);
	if(sp.ksg_w)
		whole_results = fillAndRunWholeSeriesKSG(gp, sp, allInteractions, entire_history, hAgg_history, cumulativeRotations, fullgtennIdx, posPts,
								 rng);

	if(sp.ksg_local == 1)
		localised_results = fillAndRunLocalInteractionsKSG(gp, sp, allInteractions, entire_history, hAgg_history, &avgDirections[0], flockIds, cumulativeRotations, fullgtennIdx, posPts, rng);

	if(sp.ksg_local == 2)
		localised_results = fillAndRunFlockvsParticleKSG(gp, sp, allInteractions, entire_history, hAgg_history, &avgDirections[0], flockIds, cumulativeRotations, fullgtennIdx, posPts, rng);


	/******************
	 * Write final logs
	 *****************/
	for(int i = 0; i < SIM_COUNT; ++i) {
		if(sp.sims[i] == false)
			continue;

		auto ts_stats = per_timestep_stats[i];
		auto subset_stats = subset_results[i];
		auto whole_stats = whole_results[i];
		auto localised_stats = localised_results[i];
		fmt::print(logFiles[i], "\nPer timestep based on interactions:\nI ={}\n", ts_stats.avg);
		fmt::print(logFiles[i], "\nImin ={}\n", ts_stats.min);
		fmt::print(logFiles[i], "\nImax ={}\n", ts_stats.max);
		fmt::print(logFiles[i], "\nWhole timeseries. All i,j, where i != j\nI2 ={}\n", whole_stats.avg);
		fmt::print(logFiles[i], "\nImin2 ={}\n", whole_stats.min);
		fmt::print(logFiles[i], "\nImax2 ={}\n", whole_stats.max);
		fmt::print(logFiles[i], "\nDecimation of Whole timeseries (subsets). All i,j, where i != j\nI3 ={}\n", subset_stats.avg);
		fmt::print(logFiles[i], "\nImin3 ={}\n", subset_stats.min);
		fmt::print(logFiles[i], "\nImax3 ={}\n", subset_stats.max);
		fmt::print(logFiles[i], "\nConnected flock adjusted.\nI4 ={}\n", localised_stats.avg);
		fmt::print(logFiles[i], "\nImin4 ={}\n", localised_stats.min);
		fmt::print(logFiles[i], "\nImax4 ={}\n", localised_stats.max);

		itHists.writeToFile(binFiles[i], (SIM)i);
		if(i == GTEBIN) {
			fmt::print(logFiles[i], "\n1D KL GTE\nI5 ={}\n", kl1DGTE);
		}

		if(sp.angleWrite) {
			if(fwrite(entire_history.read_h(), sizeof(double), gp.N*recordedU, binFiles[i]) != gp.N*recordedU) PEEXIT("write failed");
			if(fwrite(hAgg_history, sizeof(double), gp.N*recordedU, binFiles[i]) != gp.N*recordedU) PEEXIT("write failed");
		}
	}

	/***
	* Free memory used
	***/
	free(flockIds);
	delete[] fullgtennIdx;
	free(cumulativeRotations);
	free(rotations);
	free(hAgg_history);

	fmt::print("\nWriting file {}\n", ofile);
	endFiles(sp.sims, logFiles, binFiles, ofile.c_str());
	
	return EXIT_SUCCESS;
}

#include "../include/model.hpp"
#include "../include/stats.hpp"
#include "../include/logFileHelper.hpp"
#include "../include/version.hpp"
#include "../include/defs.hpp"
#include "../include/vpTree.hpp"
#include "../include/utils.hpp"

#include <fmt/format.h>
#include <clara.hpp>
#include <iostream>
#include <omp.h>
#include <stdlib.h>
#include <vector>
#include <random>

int sim_generate_threaded(int argc, char* argv[])
{
	std::string initialState_nonConst = "", ofile_nonConst = "incorrectName";
	int imode_nonConst = 0, threadLimit_nonConst = 1;
	GeneratorParams gp{}; // Need to set up defaults. Bools need to default to false
	{
		using namespace clara::detail;
		bool showHelp = false;
		auto cli = Help(showHelp)
			| Opt(initialState_nonConst, "filename")["-i"]["--in-file"]("Input file containing initial state. Leave blank for random initialisation.")
			| Opt(ofile_nonConst, "filename")["-o"]["--out-file"]("Base name for output data. Will prepend mibin/tebin/gtebin and append .bin or .log").required()
			| Opt(imode_nonConst, "mode")["--imode"]("initialisation mode (0 - random, 1 - random (same) or 2 - all zero)")
			| generatorParamsParser(gp);
		auto result = cli.parse(Args(argc, argv));
		if (!result) {
			fmt::print("Error in command line: {}\n", result.errorMessage());
			return EXIT_FAILURE;
		}
		if (showHelp) {
			std::cerr << cli << std::endl;
			return EXIT_SUCCESS;
		}
		gp.L = sqrt(((double)gp.N) / gp.rho);
	}
	const std::string initialState = initialState_nonConst, ofile = ofile_nonConst;
	const int imode = imode_nonConst, threadLimit = threadLimit_nonConst;

	int numThreads = omp_get_max_threads();
	omp_set_num_threads(numThreads);
	if (threadLimit > 0) {
		omp_set_num_threads(threadLimit);
		numThreads = threadLimit;
		fmt::print("Limiting threads used to {} (Using {}, max {})\n", threadLimit, numThreads, omp_get_max_threads());
	}
	fmt::print("Max threads: {} ({})\n", omp_get_max_threads(), numThreads);
	fmt::print("Using {} threads\n", numThreads);

	fmt::print("\noutput file        : {}\n", ofile);
	fmt::print("linear size (L)    : {}\n", gp.L);

	FILE* fp = fopen(ofile.c_str(), "wb");
	if(fp == NULL) PEEXIT("failed to open output file");

	gp.seed = calcSeed(gp.seed);
	std::mt19937_64 rng(gp.seed);
	fmt::print("seed : {}\n", gp.seed);

	writeBinaryHeader(fp, gp, SimParams{});

	fmt::print("Git Head: {}\n", GIT_HEAD);
	fmt::print("Git Log: {}\n", GIT_LOG);
	fmt::print("Git Diff: {}\n", GIT_DIFF);

	fmt::print("Initialising model\n");
	FlockState fs(gp.N);
	DoubleBuffer<double> sBuf(gp.N);
	double* const hAgg = (double* const)calloc(gp.N, sizeof(double));

	if(initialState.empty()) {
		initialise(imode, gp.N, gp.L, fs.h.writer(), fs.x.writer(), fs.y.writer(), sBuf.writer(), rng);
		fs.swapBuffers(); sBuf.swap();
	}
	else {
		FILE* initFile = fopen(initialState.c_str(), "rb");
		if(initFile == NULL) PEEXIT("failed to open initialisation file");

		/*auto unusedVars =*/ readBinaryHeader(initFile);
		readBinaryContents(initFile, fs);
		fclose(initFile);

		fs.swapBuffers();
		duplicate_reader_to_writer(fs.h);
		duplicate_reader_to_writer(fs.x);
		duplicate_reader_to_writer(fs.y);
	}

	//const double dt = 1;

	std::random_device r;
	std::mt19937_64 engine(r());
	std::normal_distribution<> gaussian_rng(0.0, 1.0);
	/***
	* Perform lead in phase of updates
	***/
	VpTreeEuclidean2D vicsekTree(gp.L);
	fmt::print("Performing lead-in updates\n");
	for(size_t u = 0; u < gp.S; ++u) {
		if(gp.imethod == 0) {
			if(gp.rotate_frame == 0) {
				vicsekMetricThreaded(gp.N, gp.L, gp.v, gp.eta, fs, hAgg, rng, gp.umethod, vicsekTree, 0, numThreads);
			}
		}
		/*else if(imethod == 1) {
			vicsekTopological(N, L, v, eta, dt, h, x, y, hnew, xnew, ynew, hAgg, &rng, topo_neighbours, vicsekTree, umethod);
		}*/
		else
			PEEXIT("Invalid imethod");
		
		// swap buffers
		fs.swapBuffers(); 
		sBuf.swap();

		progrep("skip", u, gp.S);
	}

	writeBinaryContents(fp, fs);

	/***
	* Free memory used
	***/
	fmt::print("\nWriting file {}\n", ofile);
	if(fclose(fp) != 0) PEEXIT("failed to close output file");

	free(hAgg);

	return EXIT_SUCCESS;
}

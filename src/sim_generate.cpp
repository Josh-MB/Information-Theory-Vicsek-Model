#include "../include/model.hpp"
#include "../include/stats.hpp"
#include "../include/logFileHelper.hpp"
#include "../include/version.hpp"
#include "../include/defs.hpp"
#include "../include/utils.hpp"

#include <fmt/format.h>
#include <clara.hpp>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <random>

int sim_generate(int argc, char* argv[])
{
	/***
	* Create variables from arguments (set to default if not in argv)
	***/
	std::string initialState_nonConst = "", ofile_nonConst = "incorrectName";
	int imode_nonConst = 0;
	GeneratorParams gp{}; // Need to set up defaults. Bools need to default to false
	{
		using namespace clara::detail;
		bool showHelp = false;
		auto cli = Help(showHelp)
			| Opt(initialState_nonConst, "filename")["-i"]["--in-file"]("Input file containing initial state. Leave blank for random initialisation.")
			| Opt(ofile_nonConst, "filename")["-o"]["--out-file"]("Output file").required()
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
		if (gp.imethod == 3) {
			gp.L = std::ceil(std::sqrt(gp.N));
			gp.N = static_cast<size_t>(gp.L * gp.L);
			fmt::print("Adjusting N to {} (Next square) and L to {} for xy model", gp.N, gp.L);
		}
	}
	const std::string initialState = initialState_nonConst, ofile = ofile_nonConst;
	const int imode = imode_nonConst;

	fmt::print("\noutput file        : {}\n", ofile);
	fmt::print("linear size (L)    : {}\n", gp.L);

	FILE* fp = fopen(ofile.c_str(), "wb");
	if (fp == NULL) PEEXIT("failed to open output file");

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
		initialise(imode, gp.N, gp.L, fs.h.writer(), fs.x.writer(), fs.y.writer(), sBuf.writer(), rng, gp.imethod == 3);
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

	const double dt = (gp.force_dt == 0 ? gp.dt_factor * 0.1 * std::sqrt(gp.J / gp.chi) : gp.force_dt);
	const double sigma = std::sqrt(2 * 2 * gp.viscosity * gp.eta);

	std::random_device r;
	std::mt19937_64 engine(r());
	std::normal_distribution<> gaussian_rng(0.0, 1.0);
	/***
	* Perform lead in phase of updates
	***/
	VpTreeEuclidean2D vicsekTree(gp.L);
	fmt::print("Performing lead-in updates\n");
	for(size_t u = 0; u < gp.S; ++u) {
		auto s = sBuf.reader();
		auto snew = sBuf.writer();

		if(gp.imethod == 0) {
			double unused_rotation;
			vicsekMetric(gp.N, gp.L, gp.v, gp.eta, fs, hAgg, rng, gp.umethod, vicsekTree, 0, 0, gp.discretise, 0, gp.rotate_frame != 0, &unused_rotation);
		}
		else if(gp.imethod == 1) {
			vicsekTopological(gp.N, gp.L, gp.v, gp.eta, dt, fs, hAgg, rng, gp.topo_neighbours, vicsekTree, gp.umethod);
		}
		else if(gp.imethod == 2) {
			double hamiltonian;
			inertialSpinModel(gp.N, gp.L, gp.v, sigma, gp.chi, gp.viscosity, gp.J, dt, gp.topo_neighbours,
							  fs, s, snew, hAgg,
							  engine, gaussian_rng, gp.umethod, vicsekTree, hamiltonian, 0, 0, 0);
		}
		else if (gp.imethod == 3) {
			xyModel(gp.N, gp.L, gp.eta, fs, hAgg, rng);
		}
		else
			PEEXIT("Invalid imethod");
		
		// swap buffers
		fs.swapBuffers(); sBuf.swap();

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

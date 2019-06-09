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

int sim_record(int argc, char* argv[])
{
	/***
	* Create variables from arguments (set to default if not in argv)
	***/
	std::string initialState_nonConst = "", ofile_nonConst = "incorrectName", dataFile_nonConst = "data.bin";
	SimParams sp{};
	{
		using namespace clara::detail;
		bool showHelp = false;
		auto cli = Help(showHelp)
			| Opt(initialState_nonConst, "filename")["-i"]["--in-file"]("Input file containing initial state. Leave blank for random initialisation.").required()
			| Opt(ofile_nonConst, "filename")["-o"]["--out-file"]("Base name for output data. Will prepend mibin/tebin/gtebin and append .bin or .log").required()
			| Opt(dataFile_nonConst, "filename")["-d"]["--data-file"]("State information for every nth step. Defaults to data.bin")
			| Opt(sp.U, "timesteps")["-U"]["--updatesteps"]("Number of update steps")
			| Opt(sp.seed, "number")["--seed"]("Random seed (0 for unpredictable)")
			| Opt(sp.reseed)["--update-initial-state"]("Overwrite the initialState file with the final state here, for use with future runs")
			| Opt(sp.record_T_steps, "gap")["--record-nth-step"]("Record every nth timestep for data analysis");
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
	const std::string initialState = initialState_nonConst, ofile = ofile_nonConst, dataFile = dataFile_nonConst;
	FILE* initFile = fopen(initialState.c_str(), "rb");
	if(initFile == NULL) PEEXIT("failed to open initialisation file");

	auto binVars = readBinaryHeader(initFile);
	GeneratorParams gp = std::get<0>(binVars);

	fmt::print("\noutput file        : {}\n", ofile);
	fmt::print("linear size (L)    : {}\n", gp.L);

	FILE* fp = fopen(ofile.c_str(), "wb");
	if(fp == NULL) PEEXIT("failed to open output file");

	sp.seed = calcSeed(sp.seed);
	std::mt19937_64 rng(sp.seed);
	fmt::print("seed : {}\n", sp.seed);

	writeBinaryHeader(fp, gp, sp);

	fmt::print("Git Head: {}\n", GIT_HEAD);
	fmt::print("Git Log: {}\n", GIT_LOG);
	fmt::print("Git Diff: {}\n", GIT_DIFF);

	fmt::print("Initialising model\n");
	FlockState fs(gp.N);
	DoubleBuffer<double> sBuf(gp.N);
	double* const hAgg = (double* const)calloc(gp.N, sizeof(double));

	readBinaryContents(initFile, fs);
	fclose(initFile);

	fs.swapBuffers();
	duplicate_reader_to_writer(fs.h);
	duplicate_reader_to_writer(fs.x);
	duplicate_reader_to_writer(fs.y);

	const double dt = gp.dt_factor * 0.1 * std::sqrt(gp.J / gp.chi);
	const double sigma = std::sqrt(2 * 2 * gp.viscosity * gp.eta);

	fmt::print("Opening data file: {}\n", dataFile);
	FILE* dataFp = fopen(dataFile.c_str(), "wb");
	if(dataFp == NULL) PEEXIT("failed to open data output file");

	if(fwrite(&gp.N, sizeof(size_t), 1, dataFp) != 1) PEEXIT("failed to write to file");
	if(fwrite(&sp.U, sizeof(size_t), 1, dataFp) != 1) PEEXIT("failed to write to file");
	if(fwrite(&gp.L, sizeof(double), 1, dataFp) != 1) PEEXIT("failed to write to file");
	if(fwrite(&sp.record_T_steps, sizeof(int), 1, dataFp) != 1) PEEXIT("failed to write to file");

	size_t timeCount = sp.U / sp.record_T_steps;
	std::vector<double> order_zx;
	std::vector<double> order_zy;
	order_zx.reserve(timeCount);
	order_zy.reserve(timeCount);
	int timeStep = 0;

	std::random_device r;
	std::mt19937_64 engine(r());
	std::normal_distribution<> gaussian_rng(0.0, 1.0);

	/***
	* Perform lead in phase of updates
	***/
	VpTreeEuclidean2D vicsekTree(gp.L);
	fmt::print("Performing lead-in updates\n");
	for(size_t u = 0; u < sp.U; ++u) {
		auto s = sBuf.reader();
		auto snew = sBuf.writer();

		if(gp.imethod == 0) {
			vicsekMetric(gp.N, gp.L, gp.v, gp.eta, fs, hAgg, rng, gp.umethod, vicsekTree, 0, 0, 0);
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
		else
			PEEXIT("Invalid imethod");
		
		if((u % sp.record_T_steps) == 0) {
			auto h = fs.h.reader(), x = fs.x.reader(), y = fs.y.reader();
			if(fwrite(&u, sizeof(size_t), 1, dataFp) != 1) PEEXIT("failed to write to file");
			double zx, zy;
			order_param(gp.N, h, &zx, &zy);
			order_zx.push_back(zx);
			order_zy.push_back(zy);
			if(fwrite(x, sizeof(double), gp.N, dataFp) != gp.N) PEEXIT("failed to write to file");
			if(fwrite(y, sizeof(double), gp.N, dataFp) != gp.N) PEEXIT("failed to write to file");
			if(fwrite(h, sizeof(double), gp.N, dataFp) != gp.N) PEEXIT("failed to write to file");
			++timeStep;
		}

		// swap buffers
		fs.swapBuffers();
		sBuf.swap();

		progrep("skip", u, sp.U);
	}
	if(fwrite(&(order_zx[0]), sizeof(double), timeCount, dataFp) != timeCount) PEEXIT("write failed");
	if(fwrite(&(order_zy[0]), sizeof(double), timeCount, dataFp) != timeCount) PEEXIT("write failed");

	fclose(dataFp);

	if(sp.reseed) {
		reseed(initialState.c_str(), gp, sp, fs);
	}

	/***
	* Free memory used
	***/
	fmt::print("\nWriting file {}\n", ofile);
	if(fclose(fp) != 0) PEEXIT("failed to close output file");
	free(hAgg);

	return EXIT_SUCCESS;
}

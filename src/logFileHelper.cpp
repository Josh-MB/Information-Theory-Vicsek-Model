#include "../include/logFileHelper.hpp"
#include "../include/version.hpp"
#include "../include/utils.hpp"
#include "../include/stats.hpp"

#include <fmt/format.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

void openFiles(bool runSims[SIM_COUNT], FILE * logFiles[SIM_COUNT], FILE * binFiles[SIM_COUNT], char const*const outputFileBase)
{
	/** Opens the log and data binary files for output opted by runSims.
	File name convention is <some/dir/structure>/<sim_name><basename><.[log|bin]>
	where the argument outputFileBase is <some/dir/structure/basename>, with directories
	being optional.
	**/
	fmt::print("Opening files\n");
	const char * pch = strrchr(outputFileBase, '/');
	//fmt::print("pch: {}, ofile: {}\n", pch, outputFileBase);
	int position = 0;
	if(pch != 0) position = (int)(pch - outputFileBase);
	fmt::print("Position: {}\n", position);
	std::string buffer;
	for(int i = 0; i < SIM_COUNT; ++i) {
		if(runSims[i] == false)	continue;
		if(position == 0)
			buffer = fmt::format("{}{}{}", simNames[i], outputFileBase, extensions[0]);
		else
			buffer = fmt::format("{:.{}}{}{}", position + 1, outputFileBase, simNames[i], pch + 1, extensions[0]);
		binFiles[i] = fopen(buffer.c_str(), "wb");
		if(binFiles[i] == NULL) PEEXIT("failed to open output file");
		if(position == 0)
			buffer = fmt::format("{}{}{}", simNames[i], outputFileBase, extensions[1]);
		else
			buffer = fmt::format("{:.{}}{}{}", position + 1, outputFileBase, simNames[i], pch + 1, extensions[1]);
		logFiles[i] = fopen(buffer.c_str(), "w");
		if(logFiles[i] == NULL) PEEXIT("failed to open output file");
	}
}

void endFiles(bool runSims[SIM_COUNT], FILE * logFiles[SIM_COUNT], FILE * binFiles[SIM_COUNT], char const*const outputFileBase)
{
	for(int i = 0; i < SIM_COUNT; i++) {
		if(runSims[i]) {
			fmt::print(logFiles[i], "\nWriting file {}{}{}\n", simNames[i], outputFileBase, extensions[0]);
			if(fclose(binFiles[i]) != 0) PEEXIT("failed to close output file");
			fclose(logFiles[i]);
		}
	}
}

void writeFileHeaders(FILE * logFiles[SIM_COUNT], FILE * binFiles[SIM_COUNT], GeneratorParams const & gp, SimParams const & sp, char const*const outputFileBase)
{
	fmt::print("Writing headers\n");
	for (int i = 0; i < SIM_COUNT; ++i) {
		if (sp.sims[i] == false)
			continue;
		fmt::print(logFiles[i], "{:<12} = {:<24} {}\n", "N", gp.N, "number of particles (density)");
		fmt::print(logFiles[i], "{:<12} = {:<24} {}\n", "rho", gp.rho, "particle interaction intensity [0.013]");
		fmt::print(logFiles[i], "{:<12} = {:<24} {}\n", "v", gp.v, "particle velocity");
		fmt::print(logFiles[i], "{:<12} = {:<24} {}\n", "eta", gp.eta, "noise");
		fmt::print(logFiles[i], "{:<12} = {:<24} {}\n", "B", sp.B, "number of bins");
		fmt::print(logFiles[i], "{:<12} = {:<24} {}\n", "U", sp.U, "number of updates");
		fmt::print(logFiles[i], "KSG dimensions for GTE: {}\n", sp.ksg_gte_dims);
		fmt::print(logFiles[i], "Hists dimensions for GTE: {}\n", sp.hist_gte_dims);
		fmt::print(logFiles[i], "{:<12} = {:<24} {}\n", "GK", GK, "number of Y neighbours/dimensions for GTE Kraskov Estimator (compile time constant)");
		fmt::print(logFiles[i], "{:<12} = {:<24} {}\n", "GKavg", sp.GKavg, "number of Y neighbours/dimensions for GTE to average over (0 = all, i.e. variable)");
		fmt::print(logFiles[i], "{:<12} = {:<24} {}\n", "topo_neighbours", gp.topo_neighbours, "number of nearest neighbours for Topological Update");
		fmt::print(logFiles[i], "{:<12} = {:<24} {}\n", "KSG_neighbours", sp.KSG_neighbours, "number of nearest neighbours for Kraskov Estimator");
		fmt::print(logFiles[i], "{:<12} = {:<24} {}\n", "ofile", outputFileBase, "output data file names. Will prepend mibin/tebin/gtebin and append .bin or .log");
		fmt::print(logFiles[i], "{:<12} = {:<24} {}\n", "sim seed used", sp.seed, "if 0 was used for seed input setting");
		fmt::print(logFiles[i], "Frame rotation: {}\n", gp.rotate_frame == 0 ? "None" : "Full");
		fmt::print(logFiles[i], "Discretised: {}\n", gp.discretise == 0 ? "False" : "True");
		fmt::print(logFiles[i], "{:<12} = {:<24} {}\n", "ksg_local", sp.ksg_local, "Connected flock setting");
		fmt::print(logFiles[i], "Using data generated using seed: {}, over {} lead-in steps\n", gp.seed, gp.S);
		fmt::print(logFiles[i], "Used {} update method\n", gp.imethod == 0 ? "metric" : "topological");
		fmt::print(logFiles[i], "{} updating\n", gp.umethod == 1 ? "Forward" : "Backward");
		if (sp.shuffle == 1)
			fmt::print(logFiles[i], "Shuffled data. TE shuffle dim: {}\n", (sp.te_shuffle_dim == 0 ? 'X' : (sp.te_shuffle_dim == 1 ? 'Y' : 'W')));
		else
			fmt::print(logFiles[i], "No shuffling\n");
		fmt::print(logFiles[i], "Recording every {}th step\n", sp.record_T_steps);
		fmt::print(logFiles[i], "{:<12} = {:<24} {}\n", "chi", gp.chi, "ISM - moment of inertia");
		fmt::print(logFiles[i], "{:<12} = {:<24} {}\n", "J", gp.J, "ISM - alignment strength");
		fmt::print(logFiles[i], "{:<12} = {:<24} {}\n", "viscosity", gp.viscosity, "ISM - viscosity coefficient");
		fmt::print(logFiles[i], "{:<12} = {:<24} {}\n", "dt_factor", gp.dt_factor, "ISM - dt scaling");
		fmt::print(logFiles[i], "{:<12} = {:<24} {}\n", "force_dt", gp.force_dt, "Force dt: 0 - use ISM parameter otherwise use force_dt, ignores dt_factor. (ISM and vicsek topo only)");

		fmt::print(logFiles[i], "Estimator : Gomez\n");

		fmt::print(logFiles[i], "Git Head: {}\n", GIT_HEAD);
		fmt::print(logFiles[i], "Git Log: {}\n", GIT_LOG);
		fmt::print(logFiles[i], "Git Diff: {}\n", GIT_DIFF);
		std::string buffer = fmt::format("vicsek {}", simNames[i]);

		magic_header(binFiles[i], buffer.c_str());

		writeBinaryHeader(binFiles[i], gp, sp);
	}
}

std::tuple<GeneratorParams, SimParams> readBinaryHeader(FILE * file)
{
	GeneratorParams gp{};
	SimParams sp{};
	if (file) {
		if (fread(&gp.N, sizeof(size_t), 1, file) != 1) PEEXIT("reading from binary failed");
		if (fread(&gp.rho, sizeof(double), 1, file) != 1) PEEXIT("reading from binary failed");
		if (fread(&gp.v, sizeof(double), 1, file) != 1) PEEXIT("reading from binary failed");
		if (fread(&gp.eta, sizeof(double), 1, file) != 1) PEEXIT("reading from binary failed");
		if (fread(&gp.L, sizeof(double), 1, file) != 1) PEEXIT("reading from binary failed");
		if (fread(&sp.B, sizeof(size_t), 1, file) != 1) PEEXIT("reading from binary failed");
		if (fread(&sp.U, sizeof(size_t), 1, file) != 1) PEEXIT("reading from binary failed");
		if (fread(&sp.ksg_gte_dims, sizeof(int), 1, file) != 1) PEEXIT("reading from binary failed");
		if (fread(&sp.hist_gte_dims, sizeof(int), 1, file) != 1) PEEXIT("reading from binary failed");
		size_t GK1;
		if (fread(&GK1, sizeof(size_t), 1, file) != 1) PEEXIT("reading from binary failed");
		if (GK1 != GK) PEEXIT("GTE K ({}) in binary file is different to the static const in this version (%zu)", GK1, GK);
		if (fread(&sp.GKavg, sizeof(size_t), 1, file) != 1) PEEXIT("reading from binary failed");
		if (fread(&gp.topo_neighbours, sizeof(size_t), 1, file) != 1) PEEXIT("reading from binary failed");
		if (fread(&sp.KSG_neighbours, sizeof(size_t), 1, file) != 1) PEEXIT("reading from binary failed");
		if (fread(&gp.seed, sizeof(uint64_t), 1, file) != 1) PEEXIT("reading from binary failed");
		if (fread(&sp.seed, sizeof(uint64_t), 1, file) != 1) PEEXIT("reading from binary failed");
		if (fread(&gp.imethod, sizeof(int), 1, file) != 1) PEEXIT("reading from binary failed");
		if (fread(&gp.umethod, sizeof(int), 1, file) != 1) PEEXIT("reading from binary failed");
		if (fread(&sp.te_shuffle_dim, sizeof(int), 1, file) != 1) PEEXIT("reading from binary failed");
		if (fread(&gp.discretise, sizeof(int), 1, file) != 1) PEEXIT("reading from binary failed");
		if (fread(&sp.ksg_local, sizeof(int), 1, file) != 1) PEEXIT("reading from binary failed");
		if (fread(&sp.record_T_steps, sizeof(int), 1, file) != 1) PEEXIT("reading from binary failed");
		if (fread(&gp.chi, sizeof(double), 1, file) != 1) PEEXIT("reading from binary failed");
		if (fread(&gp.J, sizeof(double), 1, file) != 1) PEEXIT("reading from binary failed");
		if (fread(&gp.viscosity, sizeof(double), 1, file) != 1) PEEXIT("reading from binary failed");
		if (fread(&gp.dt_factor, sizeof(double), 1, file) != 1) PEEXIT("reading from binary failed");
		if (fread(&gp.force_dt, sizeof(double), 1, file) != 1) PEEXIT("reading from binary failed");
		int tmpBool;
		if (fread(&tmpBool, sizeof(int), 1, file) != 1) PEEXIT("reading from binary failed");
		gp.rotate_frame = static_cast<bool>(tmpBool);
		if (fread(&tmpBool, sizeof(int), 1, file) != 1) PEEXIT("reading from binary failed");
		sp.shuffle = static_cast<bool>(tmpBool);
	}
	return { gp, sp };
}

void readBinaryContents(FILE * file, FlockState& fs)
{
	if(file != 0) {
		size_t N = fs.x.size();
		if(fread(fs.x.writer(), sizeof(double), N, file) != N) PEEXIT("reading x from binary failed");
		if(fread(fs.y.writer(), sizeof(double), N, file) != N) PEEXIT("reading y from binary failed");
		if(fread(fs.h.writer(), sizeof(double), N, file) != N) PEEXIT("reading h from binary failed");
	}
}

void writeBinaryHeader(FILE * file, GeneratorParams const & gp, SimParams const & sp)
{
	if (file != 0) {
		if (fwrite(&gp.N, sizeof(size_t), 1, file) != 1) PEEXIT("write failed");
		if (fwrite(&gp.rho, sizeof(double), 1, file) != 1) PEEXIT("write failed");
		if (fwrite(&gp.v, sizeof(double), 1, file) != 1) PEEXIT("write failed");
		if (fwrite(&gp.eta, sizeof(double), 1, file) != 1) PEEXIT("write failed");
		if (fwrite(&gp.L, sizeof(double), 1, file) != 1) PEEXIT("write failed");
		if (fwrite(&sp.B, sizeof(size_t), 1, file) != 1) PEEXIT("write failed");
		if (fwrite(&sp.U, sizeof(size_t), 1, file) != 1) PEEXIT("write failed");
		if (fwrite(&sp.ksg_gte_dims, sizeof(int), 1, file) != 1) PEEXIT("write failed");
		if (fwrite(&sp.hist_gte_dims, sizeof(int), 1, file) != 1) PEEXIT("write failed");
		if (fwrite(&GK, sizeof(size_t), 1, file) != 1) PEEXIT("write failed");
		if (fwrite(&sp.GKavg, sizeof(size_t), 1, file) != 1) PEEXIT("write failed");
		if (fwrite(&gp.topo_neighbours, sizeof(size_t), 1, file) != 1) PEEXIT("write failed");
		if (fwrite(&sp.KSG_neighbours, sizeof(size_t), 1, file) != 1) PEEXIT("write failed");
		if (fwrite(&gp.seed, sizeof(uint64_t), 1, file) != 1) PEEXIT("write failed");
		if (fwrite(&sp.seed, sizeof(uint64_t), 1, file) != 1) PEEXIT("write failed");
		if (fwrite(&gp.imethod, sizeof(int), 1, file) != 1) PEEXIT("write failed");
		if (fwrite(&gp.umethod, sizeof(int), 1, file) != 1) PEEXIT("write failed");
		if (fwrite(&sp.te_shuffle_dim, sizeof(int), 1, file) != 1) PEEXIT("write failed");
		if (fwrite(&gp.discretise, sizeof(int), 1, file) != 1) PEEXIT("write failed");
		if (fwrite(&sp.ksg_local, sizeof(int), 1, file) != 1) PEEXIT("write failed");
		if (fwrite(&sp.record_T_steps, sizeof(int), 1, file) != 1) PEEXIT("write failed");
		if (fwrite(&gp.chi, sizeof(double), 1, file) != 1) PEEXIT("write failed");
		if (fwrite(&gp.J, sizeof(double), 1, file) != 1) PEEXIT("write failed");
		if (fwrite(&gp.viscosity, sizeof(double), 1, file) != 1) PEEXIT("write failed");
		if (fwrite(&gp.dt_factor, sizeof(double), 1, file) != 1) PEEXIT("write failed");
		if (fwrite(&gp.force_dt, sizeof(double), 1, file) != 1) PEEXIT("write failed");
		int tmpBool = static_cast<int>(gp.rotate_frame);
		if (fwrite(&tmpBool, sizeof(int), 1, file) != 1) PEEXIT("write failed");
		tmpBool = static_cast<int>(sp.shuffle);
		if (fwrite(&tmpBool, sizeof(int), 1, file) != 1) PEEXIT("write failed");
	}
}

void writeBinaryContents(FILE * file, FlockState const& fs)
{
	if(file != 0) {
		size_t N = fs.x.size();
		if(fwrite(fs.x.reader(), sizeof(double), N, file) != N) PEEXIT("writing x to file failed");
		if(fwrite(fs.y.reader(), sizeof(double), N, file) != N) PEEXIT("reading y to file failed");
		if(fwrite(fs.h.reader(), sizeof(double), N, file) != N) PEEXIT("reading h to file failed");
	}
}

void calc_and_write_order_param(size_t N, double const * const h, double rot, bool sims[SIM_COUNT], FILE* binFiles[SIM_COUNT])
{
	double zx, zy;
	order_param(N, h, &zx, &zy, rot);

	for (int i = 0; i < SIM_COUNT; ++i) {
		if (sims[i] == false)
			continue;
		if (fwrite(&zx, sizeof(double), 1, binFiles[i]) != 1) PEEXIT("write failed");
		if (fwrite(&zy, sizeof(double), 1, binFiles[i]) != 1) PEEXIT("write failed");
	}
}

void reseed(char const * const initialState, GeneratorParams const & gp, SimParams const & sp, FlockState const & fs)
{
	FILE* fp = fopen(initialState, "wb");
	fmt::print("Reseeding\n");
	if (fp == NULL) PEEXIT("failed to open initialstate file for reseeding");
	writeBinaryHeader(fp, gp, sp);
	writeBinaryContents(fp, fs);
	if (fclose(fp) != 0) PEEXIT("failed to close output file");
}

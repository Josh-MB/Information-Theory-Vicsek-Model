#include "../include/anim_utils.hpp"
#include "../include/defs.hpp"
#include "../include/utils.hpp"

#include <fmt/format.h>
#include <clara.hpp>
#include <iostream>
#include <list>
#include <stdlib.h>

#ifdef ___OSX___
	#define M_PI 3.14159265358979323846
#endif

int sim_visualize(int argc, char* argv[])
{
	std::string gnuplot_args = "''", outputDir_nonConst = "frames", dataFile_nonConst = "data.bin";
	int filetype_nonConst = 0;
	{
		using namespace clara::detail;
		bool showHelp = false;
		auto cli = Help(showHelp)
			| Opt(gnuplot_args, "'<args>'")["--gnuplot-args"]("Arguments for gnuplot command")
			| Opt(filetype_nonConst, "filetype")["--filetype"]("0 = png, 1 = eps")
			| Opt(outputDir_nonConst, "directory")["--out-dir"]("Output directory")
			| Opt(dataFile_nonConst, "filename")["-i"]["--in-file"]["--data-file"]("Input file");
		auto result = cli.parse(Args(argc, argv));
		if (!result) {
			fmt::print("Error in command line: {}\n", result.errorMessage());
			return EXIT_FAILURE;
		}
		if (!gnuplot_args.empty() && (gnuplot_args.front() != '\'' || gnuplot_args.back() != '\'')) {
			fmt::print("Error in command line: --gnuplot-args must be single-quoted string");
			return EXIT_FAILURE;
		}
		if (showHelp) {
			std::cerr << cli << std::endl;
			return EXIT_SUCCESS;
		}
	}

	const std::string gp_args = gnuplot_args.empty() ? "" : gnuplot_args.substr(1, gnuplot_args.size() - 2);
	const std::string gpcmd = "gnuplot " + gp_args;
	const std::string outputDir = outputDir_nonConst, dataFile = dataFile_nonConst;
	const int filetype = filetype_nonConst;

	make_dir(outputDir);
	FILE *fp = fopen(dataFile.c_str(), "rb");
	if(fp == NULL) PEEXIT("failed to open data input file");

	size_t N, U;
	double L;
	int record_T_steps = 1000;
	if(fread(&N, sizeof(size_t), 1, fp) != 1) PEEXIT("failed to read N");
	if(fread(&U, sizeof(size_t), 1, fp) != 1) PEEXIT("failed to read U");
	if(fread(&L, sizeof(double), 1, fp) != 1) PEEXIT("failed to read U");
	if(fread(&record_T_steps, sizeof(int), 1, fp) != 1) PEEXIT("failed to read U");
	const double arrowlen = 0.03*L;

	double* const x = (double* const)malloc(N * sizeof(double));
	double* const y = (double* const)malloc(N * sizeof(double));
	double* const h = (double* const)malloc(N * sizeof(double));

	int* const r = (int* const)calloc(N, sizeof(int));
	int* const g = (int* const)calloc(N, sizeof(int));
	int* const b = (int* const)calloc(N, sizeof(int));

	FILE* gpipe = popen(gpcmd.c_str(), "w");
	if(gpipe == NULL)                     PEEXIT("failed to open pipe to '{}'\n", gpcmd);
	if(setvbuf(gpipe, NULL, _IOLBF, 0) != 0) PEEXIT("failed to line-buffer pipe to '{}'\n", gpcmd);

	size_t timeCount = U / record_T_steps;
	
	for(size_t u = 0; u < timeCount; ++u) {
		size_t t = u;
		if(fread(&t, sizeof(size_t), 1, fp) != 1) PEEXIT("failed to read t");
		if(fread(x, sizeof(double), N, fp) != N) PEEXIT("failed to read x");
		if(fread(y, sizeof(double), N, fp) != N) PEEXIT("failed to read y");
		if(fread(h, sizeof(double), N, fp) != N) PEEXIT("failed to read h");

		//Set particle colours based on angle
		for(size_t i = 0; i < N; ++i)
		{
			double val = (h[i] + PI) / (PI);
			hslToRgb(val - (int)val, 0.75, 0.5, r[i], g[i], b[i]);
		}

		//Set up gnuplot
		std::string frameFileName;
		if(filetype == 1) {
			fmt::print(gpipe, "set terminal postscript eps enhanced size 6in,6in\n");
			frameFileName = fmt::format("{}/frame{:04}.ps", outputDir, u);
		}
		else if(filetype == 0) {
			fmt::print(gpipe, "set term png size 1600,1000\n");
			frameFileName = fmt::format("{}/frame{:04}.png", outputDir, u);
		}

		fmt::print(gpipe, "set out \"{}\"\n", frameFileName);
		fmt::print(gpipe, "set termopt enhanced\n");

		//Prepare plot
		fmt::print(gpipe, "set title 'Vicsek Visualisation: N={}, t={}'\n", N, u);
		fmt::print(gpipe, "unset key\n");
		fmt::print(gpipe, "set xr [0:{}]\n", L);
		fmt::print(gpipe, "set yr [0:{}]\n", L);
		fmt::print(gpipe, "set size square\n");
		fmt::print(gpipe, "unset auto\n");
		fmt::print(gpipe, "unset xlabel\n");
		fmt::print(gpipe, "unset ylabel\n");
		fmt::print(gpipe, "rgb(r,g,b) = int(r)*65536 + int(g)*256 + int(b)\n");
		
		//Plot data
		fmt::print(gpipe, "plot '-' u ($2):($3):({0}*cos($1)):({0}*sin($1)):(rgb($4,$5,$6)) w vectors notitle lc rgb variable\n", arrowlen);
		for(size_t i = 0; i<N; ++i) { // for each particle
			fmt::print(gpipe, "{:20.16} {:20.16} {:20.16} {} {} {}\n", h[i], x[i], y[i], r[i], g[i], b[i]);
		}
		fmt::print(gpipe, "e\n");
	}
	std::vector<double> order_zx(timeCount);
	std::vector<double> order_zy(timeCount);
	if(fread(&(order_zx[0]), sizeof(double), timeCount, fp) != timeCount) PEEXIT("failed to read zx");
	if(fread(&(order_zy[0]), sizeof(double), timeCount, fp) != timeCount) PEEXIT("failed to read zy");

	free(r);
	free(g);
	free(b);

	free(x);
	free(y);
	free(h);

	fclose(fp);
	fclose(gpipe);

	return EXIT_SUCCESS;
}

#include "../include/model.hpp"
#include "../include/stats.hpp"
#include "../include/vpUtil.hpp"
#include "../include/dataSet.hpp"
#include "../include/metrics.hpp"
#include "../include/metricTests.hpp"
#include "../include/connectedFlocks.hpp"
#include "../include/anim_utils.hpp"
#include "../include/utils.hpp"

#include <clara.hpp>
#include <iostream>
#include <list>
#include <random>
#include <stdlib.h>

#if defined(___OSX___) || defined(WIN32)
#define M_PI 3.14159265358979323846
#endif

int sim_anim(int argc, char* argv[])
{
	std::string gnuplot_args = "''", outputDir_nonConst = "frames";
	int imode_nonConst = 0, renderFrames_nonConst = 10, outputNum_nonConst = 0;
	int filetype_nonConst = 0, orderWindowSize_nonConst = 1;
	GeneratorParams gp{};
	SimParams sp{};
	{
		using namespace clara::detail;
		bool showHelp = false;
		auto cli = Help(showHelp)
			| Opt(imode_nonConst, "mode")["--imode"]("initialisation mode (0 - random, 1 - random (same) or 2 - all zero)")
			| Opt(gnuplot_args, "'<args>'")["--gnuplot-args"]("Arguments for gnuplot command")
			| Opt(renderFrames_nonConst, "gap")["--render-nth-step"]("Render every nth timestep")
			| Opt(outputNum_nonConst, "id")["--output-num"]("Output number for order param data")
			| Opt(filetype_nonConst, "filetype")["--filetype"]("0 = png, 1 = eps")
			| Opt(outputDir_nonConst, "directory")["--out-dir"]("Output directory")
			| Opt(orderWindowSize_nonConst, "size")["--order-window-size"]("Time window size to track order history")
			| generatorParamsParser(gp)
			| simParamsParser(sp);
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
	const int imode = imode_nonConst, renderFrames = renderFrames_nonConst, outputNum = outputNum_nonConst;
	const int filetype = filetype_nonConst, orderWindowSize = orderWindowSize_nonConst;
	const std::string gp_args = gnuplot_args.empty() ? "" : gnuplot_args.substr(1, gnuplot_args.size() - 2);
	const std::string gpcmd = "gnuplot " + gp_args;
	const std::string outputDir = outputDir_nonConst;

	make_dir(outputDir);

	gp.L = sqrt(((double)gp.N) / gp.rho);
	if (gp.imethod == 3) {
		gp.L = std::ceil(std::sqrt(gp.N));
		gp.N = static_cast<size_t>(gp.L * gp.L);
		fmt::print("Adjusting N to {} (Next square) and L to {} for xy model", gp.N, gp.L);
	}
	const double n = M_PI*gp.rho;
	fmt::print("\nL = {}\nn = {}\n\n",gp.L,n);

	gp.seed = calcSeed(gp.seed);
	std::mt19937_64 rng(gp.seed);
	fmt::print("seed : {}\n", gp.seed);

	std::mt19937_64 topoRng(gp.seed);

	srand(1);

	const double arrowlen = 0.03*gp.L;
	const int d = 2;
	
	const double dt = gp.dt_factor * 0.1 * std::sqrt(gp.J / gp.chi);
	
	FlockState fs(gp.N);
	DoubleBuffer<double> sBuf(gp.N);
	double* const hAgg = (double* const)calloc(gp.N,sizeof(double));

	const int orderN = 1000;
	int orderNCurr = 0;
	double* const zx = (double* const)calloc(orderN,sizeof(double));
	double* const zy = (double* const)calloc(orderN,sizeof(double));
	double* const order = (double* const)calloc(orderN, sizeof(double));
	int* const rc = (int* const)calloc(gp.N, sizeof(int));
	int* const gc = (int* const)calloc(gp.N, sizeof(int));
	int* const bc = (int* const)calloc(gp.N, sizeof(int));
	
	double* const zxWindowBuffer = (double* const)calloc(orderWindowSize, sizeof(double));
	double* const zyWindowBuffer = (double* const)calloc(orderWindowSize, sizeof(double));
	double* const zxWindowRender = (double* const)calloc(orderWindowSize, sizeof(double));
	double* const zyWindowRender = (double* const)calloc(orderWindowSize, sizeof(double));
	int currentzWindowIdx = -1;
	int* const orderrc = (int* const)calloc(orderWindowSize,sizeof(int));
	int* const ordergc = (int* const)calloc(orderWindowSize,sizeof(int));
	int* const orderbc = (int* const)calloc(orderWindowSize,sizeof(int));
	
	for(int i = 0; i < orderWindowSize; i++)
	{
		double t = (double)i/orderWindowSize;
		hslToRgb(DEGTORAD((1-t)*240 + t*360)/TWOPI, 0.75, 0.5, orderrc[i], ordergc[i], orderbc[i]);
	}

	initialise(imode, gp.N, gp.L,fs.h.writer(),fs.x.writer(),fs.y.writer(),sBuf.writer(), rng, gp.imethod == 3);
	fs.swapBuffers();
	sBuf.swap();
	
	std::string frameFileName;
	double eta1 = 2.0;
	//double hU = (double)U / 2;
	FILE* gpipe = popen(gpcmd.c_str(), "w");
	if(gpipe == NULL)                     PEEXIT("failed to open pipe to '{}'\n", gpcmd);
	if(setvbuf(gpipe, NULL, _IOLBF, 0) != 0) PEEXIT("failed to line-buffer pipe to '{}'\n", gpcmd);

	size_t frameCount = 0;
	VpTreeEuclidean2D vicsekTree(gp.L, &topoRng);
	std::vector<int> allInteractions;
	FILE* fpOrder;
	std::string buffer = fmt::format("{}/orderParam{}.bin", outputDir, outputNum);
	fpOrder = fopen(buffer.c_str(), "wb");

	FILE* fpFlockSizes;
	buffer = fmt::format("{}/flockSize{}.bin", outputDir, outputNum);
	fpFlockSizes = fopen(buffer.c_str(), "wb");
	
	int lrun = 1000;
	int mrun = 20000;
	int hrun = 50000;
						// 0,	 1,    2,    3,    4,    5,    6,    7,    8,    9,   10,   11,   12,   13,   14,   15,   16,   17,   18,   19,   20,   21,   22,   23,   24,   25,   26
	double etaVar[] = { 6.28, 6.00, 5.50, 5.00, 4.50, 4.00, 3.50, 3.00, 2.50, 2.00, 1.70, 1.60, 1.50, 1.40, 1.30, 1.20, 1.10, 1.00, 0.90, 0.80, 0.70, 0.60, 0.50, 0.40, 0.30, 0.20, 0.10 };
	//double etaVar[] = { 1000000000000, 100000000000, 10000000000, 1000000000, 100000000, 10000000, 1000000, 100000, 10000, 1000, 100, 10, 1, .01, .001, .0001 };
	int etalen = sizeof(etaVar)/sizeof(double);
	fmt::print("Size of etaVar: {}\n", etalen);
	//int steps[] = {lrun, lrun, lrun, lrun, lrun, mrun, mrun, mrun, mrun, mrun, mrun, hrun, hrun, hrun, hrun, hrun };
	//int steps[] = {     lrun, lrun, lrun, lrun, lrun, lrun, lrun, lrun, lrun, lrun, mrun, mrun,	mrun, mrun, mrun, mrun, mrun, mrun, hrun, hrun, hrun, hrun, hrun, hrun, hrun, hrun, hrun };
	double MIDiag = 0;
	if(fwrite(&sp.U, sizeof(size_t), 1, fpOrder) != 1) PEEXIT("write failed");
	if(fwrite(&sp.U, sizeof(size_t), 1, fpFlockSizes) != 1) PEEXIT("write failed");
	size_t trueT = 0;
	//int miCount = 0;

	std::vector<double> rVals;
	const double maxR = 0.5 * sqrt(2 * gp.L*gp.L);
	for(double r = 0.0; r < maxR; r += 1.0) {
		rVals.push_back(r);
	}

	ConnectedSets connectedSet(gp.N);
	std::vector<double> avgDirections;

	double* const cosVelocities = (double*const)calloc(gp.N*sp.U, sizeof(double));
	double* const accelerations = (double*const)calloc(gp.N*sp.U, sizeof(double));

	std::vector<double> orderPlot, susceptPlot;

	std::vector<double> tmpOrder;

	std::random_device r;
	std::mt19937_64 engine(r());
	std::normal_distribution<> gaussian_rng(0.0, 1.0);
	//bool startDrawing = false;
	//for(size_t run = 0; run < 27; ++run) {
	//for(size_t run = 0; run < etalen; ++run) {
	{
		auto h = fs.h.reader(), x = fs.x.reader(), y = fs.y.reader(), s = sBuf.reader();
		auto hnew = fs.h.writer(), snew = sBuf.writer();
		//eta1 = etaVar[run];
		eta1 = gp.eta;
		tmpOrder.clear();
		const double sigma = std::sqrt(2 * d * gp.viscosity * eta1);
		size_t skipSteps = mrun;
		if(gp.eta >= 3.00) skipSteps = lrun;
		else if(gp.eta <= 1.00) skipSteps = hrun;
		//for(size_t u = 0; u < steps[run] + U; ++u, ++trueT) { // for each update
		for(size_t u = 0; u < skipSteps + sp.U; ++u, ++trueT) { // for each update
			allInteractions.clear();
			connectedSet.clear();
			
			double rotation;

			//Run model
			if(gp.imethod == 0) {
				vicsekMetric(gp.N, gp.L, gp.v, eta1, fs, hAgg, rng, gp.umethod, vicsekTree, &(allInteractions), 0, gp.discretise, 0, gp.rotate_frame != 0, &rotation);
			}
			else if(gp.imethod == 1) {
				vicsekTopological(gp.N, gp.L, gp.v, eta1, dt,fs, hAgg, rng, gp.topo_neighbours, vicsekTree, gp.umethod, &(allInteractions), 0);
			}
			else if(gp.imethod == 2) {
				//extendedVicsekMetric(N, L, v, eta1, h, x, y, hnew, xnew, ynew, hAgg, &rng, umethod, cohesionStrength, separationStrength, vicsekTree, &(allInteractions), 0, discretise, &connectedSet);
				double hamiltonian;
				inertialSpinModel(gp.N, gp.L, gp.v, sigma, gp.chi, gp.viscosity, gp.J, dt, gp.topo_neighbours, fs, s, snew, hAgg, engine, gaussian_rng, gp.umethod, vicsekTree, hamiltonian, &(allInteractions), 0, 0);
			}
			else if (gp.imethod == 3) {
				xyModel(gp.N, gp.L, gp.eta, fs, hAgg, rng);
			}
			else
				PEEXIT("Invalid imethod");

			if(gp.imethod == 1 && gp.N*gp.topo_neighbours != allInteractions.size() / 2)
				fmt::print("Something wrong with interactions size: {} vs {}\n", gp.N*gp.topo_neighbours, (int)allInteractions.size());
			
			//Set particle colours based on angle
			for(size_t i = 0; i < gp.N; ++i)
			{
				//int r, g, b;
				double val = (h[i] + PI) / (PI);
				hslToRgb(val - (int)val, 0.75, 0.5, rc[i], gc[i], bc[i]);
			}

			//Calculate order and susceptibility
			double izx, izy;
			order_param(gp.N, hnew, &izx, &izy);
			double instOrder = sqrt(izx*izx + izy*izy);
			if(fwrite(&izx, sizeof(double), 1, fpOrder) != 1) PEEXIT("write failed");
			if(fwrite(&izy, sizeof(double), 1, fpOrder) != 1) PEEXIT("write failed");
			order[u % orderN] = instOrder;
			orderNCurr = std::min(orderN, (int)u);
			double avgOrder = 0.0;
			for(int i = 0; i < orderNCurr; ++i) {
				avgOrder += order[i];
			}
			avgOrder /= orderNCurr;
			double susceptibility = 0.0;
			for(int i = 0; i < orderNCurr; ++i) {
				double a = (order[i] - avgOrder);
				susceptibility += (a*a);
			}
			susceptibility /= orderNCurr;
			
			
			zxWindowBuffer[u % orderWindowSize] = izx;
			zyWindowBuffer[u % orderWindowSize] = izy;
			currentzWindowIdx = (currentzWindowIdx + 1) % orderWindowSize;
			if(currentzWindowIdx != (int)u % orderWindowSize) PEEXIT("Messed up window idx");
			
			//if(u >= steps[run] + startU && (u % renderFrames) == 0) {
			if(u >= skipSteps && (u % renderFrames) == 0) {
				tmpOrder.emplace_back(instOrder);
				if(filetype == 1) {
					fmt::print(gpipe, "set terminal postscript eps enhanced size 6in,6in\n");
					frameFileName = fmt::format("{}/frame{:04}.ps", outputDir, frameCount);
				}
				else if(filetype == 0) {
					//fmt::print(gpipe, "set term png size 1280,750\n");
					fmt::print(gpipe, "set term png size 1600,1000\n");
					frameFileName = fmt::format("{}/frame{:04}.png", outputDir, frameCount);
				}
				frameCount++;
				fmt::print(gpipe, "set out \"{}\"\n", frameFileName);
				fmt::print(gpipe, "set termopt enhanced\n");
				fmt::print(gpipe, "set multiplot layout 1,2 title \"{}, N = {}, rho = {} (L = {}), v = {}, T_K={}, eta = {:.2}, MI = {:.2} Bits, t= {}, Order(t)={:.2}, Order(avg {})={:.2}, Susceptibility={:.2}\\nJ={:.2}, chi={:.2}, eta(viscosity)={:.2}, dt={}\"\n",
					(gp.imethod == 1 ? "Topological" : "Metric"), gp.N, gp.rho, gp.L, gp.v, gp.topo_neighbours, eta1, MIDiag, trueT, instOrder, orderNCurr, avgOrder, susceptibility, gp.J, gp.chi, gp.viscosity, dt);
				fmt::print(gpipe, "set lmargin 8\n");
				fmt::print(gpipe, "set rmargin 0\n");
				fmt::print(gpipe, "set tmargin 0\n");
				fmt::print(gpipe, "set bmargin 0\n");

				fmt::print(gpipe, "set title 'Order param history (t={})'\n", orderWindowSize);
				fmt::print(gpipe, "unset key\n");
				fmt::print(gpipe, "set xr [-1:1]\n");
				fmt::print(gpipe, "set yr [-1:1]\n");
				fmt::print(gpipe, "set size square\n");
				fmt::print(gpipe, "unset auto\n");
				fmt::print(gpipe, "unset border\n");
				fmt::print(gpipe, "unset tics\n");
				fmt::print(gpipe, "set object circle at 0,0 size 1\n");
				fmt::print(gpipe, "rgb(r,g,b) = int(r)*65536 + int(g)*256 + int(b)\n");
				fmt::print(gpipe, "plot '-' u ($1):($2):(rgb($3,$4,$5)) w lines notitle lc rgb variable\n");

				//Prep render buffer
				int numToCopy = (orderWindowSize - currentzWindowIdx - 1);
				if(currentzWindowIdx < (orderWindowSize - 1)) {
					memcpy(&(zxWindowRender[0]), &(zxWindowBuffer[currentzWindowIdx + 1]), numToCopy * sizeof(double));
					memcpy(&(zyWindowRender[0]), &(zyWindowBuffer[currentzWindowIdx + 1]), numToCopy * sizeof(double));
				}
				memcpy(&(zxWindowRender[numToCopy]), &(zxWindowBuffer[0]), (currentzWindowIdx + 1) * sizeof(double));
				memcpy(&(zyWindowRender[numToCopy]), &(zyWindowBuffer[0]), (currentzWindowIdx + 1) * sizeof(double));
				//fmt::print("{} {}\n", numToCopy, currentzWindowIdx);
				//zxWindowRender[orderWindowSize - 1] = -0.5;
				//0 1 2 3 4 5 6 7 8 9
				//^  
				//num= 10-0-1=9
				//copy 9 starting at [1] into [0]
				//copy 1 starting from [0] into [num+1]
				for(int i = 0; i < orderWindowSize; ++i) {
					fmt::print(gpipe, "{:20.16} {:20.16} {} {} {}\n", zxWindowRender[i], zyWindowRender[i], orderrc[i], ordergc[i], orderbc[i]);
				}

				fmt::print(gpipe, "e\n");

				fmt::print(gpipe, "set title 'Vicsek Model'\n");

				fmt::print(gpipe, "set size square\n");
				fmt::print(gpipe, "set xr [0:{:.16}]\n", gp.L);
				fmt::print(gpipe, "set yr [0:{:.16}]\n", gp.L);
				fmt::print(gpipe, "unset xlabel\n");
				fmt::print(gpipe, "unset ylabel\n");
				fmt::print(gpipe, "unset object\n");
				fmt::print(gpipe, "set border\n");
				fmt::print(gpipe, "set tics\n");
				fmt::print(gpipe, "rgb(r,g,b) = int(r)*65536 + int(g)*256 + int(b)\n");
				fmt::print(gpipe, "plot '-' u ($2):($3):({0}*cos($1)):({0}*sin($1)):(rgb($4,$5,$6)) w vectors notitle lc rgb variable\n", arrowlen);
				fprintv(gpipe, gp.N, h, x, y, rc, gc, bc);
				fmt::print(gpipe, "e\n");

				fmt::print(gpipe, "unset multiplot\n");
				fmt::print(gpipe, "unset out\n");
			}

			// swap buffers
			fs.swapBuffers();
			sBuf.swap();
		}

		double avgOrder = 0.0;
		for(size_t i = 0; i < tmpOrder.size(); ++i) {
			avgOrder += tmpOrder[i];
		}
		avgOrder /= (double)tmpOrder.size();
		double susceptibility = 0.0;
		for(size_t i = 0; i < tmpOrder.size(); ++i) {
			double a = (tmpOrder[i] - avgOrder);
			susceptibility += (a*a);
		}
		susceptibility /= (double)tmpOrder.size();

		susceptPlot.emplace_back(susceptibility);
		orderPlot.emplace_back(avgOrder);
	}

	fmt::print(gpipe, "set term png size 1600,1000\n");
	frameFileName = fmt::format("{}/order_and_suscept.png", outputDir);
	fmt::print(gpipe, "set out \"{}\"\n", frameFileName);
	fmt::print(gpipe, "set termopt enhanced\n");
	fmt::print(gpipe, "set multiplot layout 1,2 title \"Order and Susceptibility\"\n");

	fmt::print(gpipe, "set lmargin 8\n");
	fmt::print(gpipe, "set rmargin 0\n");
	fmt::print(gpipe, "set tmargin 0\n");
	fmt::print(gpipe, "set bmargin 0\n");

	double minx = etaVar[etalen - 1] * 0.1, maxx = etaVar[0] * 10;
	fmt::print(gpipe, "set title 'Average Order'\n", orderWindowSize);
	fmt::print(gpipe, "unset key\n");
	fmt::print(gpipe, "set xr [{}:{}]\n", minx, maxx);
	fmt::print(gpipe, "set yr [0:1]\n");
	fmt::print(gpipe, "set size square\n");
	fmt::print(gpipe, "unset auto\n");
	fmt::print(gpipe, "set logscale x 10\n");
	fmt::print(gpipe, "set format x '%.0tE%+T'\n");

	fmt::print(gpipe, "plot '-' u ($1):($2) w lines notitle\n");
	for(unsigned int i = 0; i < orderPlot.size(); ++i) {
		fmt::print(gpipe, "{:20.16} {:20.16}\n", etaVar[i], orderPlot[i]);
	}
	fmt::print(gpipe, "e\n");

	double maxSuscept = *std::max_element(susceptPlot.begin(), susceptPlot.end());
	fmt::print(gpipe, "set title 'Susceptibility'\n", orderWindowSize);
	fmt::print(gpipe, "unset key\n");
	fmt::print(gpipe, "set xr [{}:{}]\n", minx, maxx);
	fmt::print(gpipe, "set yr [0:{}]\n", maxSuscept);
	fmt::print(gpipe, "set size square\n");
	fmt::print(gpipe, "unset auto\n");
	fmt::print(gpipe, "set logscale x 10\n");
	fmt::print(gpipe, "set format x '%.0tE%+T'\n");

	fmt::print(gpipe, "plot '-' u ($1):($2) w lines notitle\n");
	for(unsigned int i = 0; i < susceptPlot.size(); ++i) {
		fmt::print(gpipe, "{:20.16} {:20.16}\n", etaVar[i], susceptPlot[i]);
	}
	fmt::print(gpipe, "e\n");
	fmt::print(gpipe, "unset multiplot\n");
	fmt::print(gpipe, "unset out\n");

	pclose(gpipe);
	fclose(fpOrder);
	fclose(fpFlockSizes);

	/*if (pclose(gpipe1) == -1) PEEXIT("failed to close pipe 1 to '{}'\n",gpcmd);
	if (pclose(gpipe2) == -1) PEEXIT("failed to close pipe 2 to '{}'\n",gpcmd);
	if (pclose(gpipe3) == -1) PEEXIT("failed to close pipe 3 to '{}'\n",gpcmd);*/
	free(cosVelocities);
	free(accelerations);

	free(zxWindowRender);
	free(zyWindowRender);
	free(zxWindowBuffer);
	free(zyWindowBuffer);
	free(orderrc);
	free(ordergc);
	free(orderbc);
	
	free(hAgg);

	free(zx);
	free(zy);
	free(order);

	return EXIT_SUCCESS;
}

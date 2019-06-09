#include "../include/model.hpp"
#include "../include/stats.hpp"
#include "../include/vpUtil.hpp"
#include "../include/dataSet.hpp"
#include "../include/metrics.hpp"
#include "../include/metricTests.hpp"
#include "../include/connectedFlocks.hpp"
#include "../include/version.hpp"
#include "../include/anim_utils.hpp"

#include <fmt/format.h>
#include <list>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <stdlib.h>

#ifdef ___OSX___
#define M_PI 3.14159265358979323846
#endif

int sim_singleflock(int argc, char* argv[])
{
	size_t N = 10000, target_flock_size = 1000, U_form = 500000, U_break = 500000, topo_neighbours = 6,
		KSG_neighbours = 3;
	double break_apart_threshold = 0.7, eta_threshold = 3.0, rho = 1.0, v = 0.1, dt = 1.;
	int start_eta = 0, imode = 0, imethod = 0, umethod = 0, windowSize = 50, render_formation_frames = 1,
		render_break_frames = 1;
	std::string eta_file = "eta_steps.txt", gnuplot_args = "''", save_state = "", load_state = "";
	bool rotate_frame = false, calc_order_just_lcc = false, calc_cooling_gte = false, dont_cull = false,
		calc_metric_allparticles = false, calc_mi_against_flock = false, calc_gte = false;
	uint64_t seed = 0;
	{
		using namespace clara::detail;
		bool showHelp = false;
		auto cli = Help(showHelp)
			| Opt(gnuplot_args, "'<args>'")["--gnuplot-args"]("Arguments for gnuplot command")
			| Opt(eta_file, "filename")["--eta-file"]("File with list of noise steps. Default = eta_steps.txt")
			| Opt(save_state, "filename")["-o"]["--out-file"]("When flock is formed, save the data here")
			| Opt(load_state, "filename")["-i"]["--in-file"]("Load already formed flock")
			| Opt(N, "num")["-N"]["--num-particles"]("Number of particles")
			| Opt(target_flock_size, "num")["--target-flock-size"]("Target size for a connected flock")
			| Opt(U_form, "timesteps")["--time-to-form"]("Number of timesteps to run to allow flock to form")
			| Opt(U_break, "timesteps")["--time-to-break"]("Number of timesteps to run to watch flock break apart")
			| Opt(topo_neighbours, "num")["--topo-neighbours"]("Number of neighbours for topological updating")
			| Opt(KSG_neighbours, "number")["--KSG-neighbours"]("How many neighbours to use for the KSG estimator")
			| Opt(break_apart_threshold, "percentage")["--disperse-threshold"]("Threshold for declaring a single flock has broken apart (based on --target-flock-size)")
			| Opt(eta_threshold, "percentage")["--eta-threshold"]("Flocks formed above this eta are ignored")
			| Opt(rho, "density")["--rho"]("Particle interaction intensity")
			| Opt(v, "velocity")["-v"]["--velocity"]("Particle velocity")
			| Opt(dt, "factor")["--dt"]("Timestep delta")
			| Opt(start_eta, "index")["--start-eta-index"]("Start partway through --eta-file")
			| Opt(imode, "mode")["--imode"]("initialisation mode (0 - random, 1 - random (same) or 2 - all zero)")
			| Opt(imethod, "method")["--imethod"]("Interaction Method, 0 = metric, 1 = topological")
			| Opt(umethod, "method")["--umethod"]("Update Method, 0 = backwards, 1 = forwards")
			| Opt(windowSize, "timesteps")["--window-size"]("Window size for calculating metrics")
			| Opt(render_formation_frames, "gap")["--render-nth-step-form"]("Render every nth step during formation stage. 0 to disable")
			| Opt(render_break_frames, "gap")["--render-nth-step-break"]("Render every nth step during break apart stage. 0 to disable")
			| Opt(rotate_frame)["--rotate-frame"]("Rotate the reference frame between updates")
			| Opt(calc_order_just_lcc)["--calc-order-just-largest-flock"]("Measure order over only the largest connected flock")
			| Opt(calc_cooling_gte)["--calc-cooling-gte"]("Measure GTE of simulations during cooling stage. Sanity check mostly")
			| Opt(dont_cull)["--simulate-with-all"]("Whether or not to cull all points outside the largest connected flock after the cooling stage")
			| Opt(calc_metric_allparticles)["--calc-metric-with-all"]("Whether or not to consider all points when calculating metrics, otherwise just uses the largest connected flock")
			| Opt(calc_mi_against_flock)["--calc-mi-against-flock"]("Measure MI between headings of particles and heading of largest connected flock")
			| Opt(calc_gte)["--calc-gte"]("Measure GTE of flock")
			| Opt(rotate_frame)["--rotate-frame"]("Rotate the reference frame between updates")
			| Opt(seed, "number")["--seed"]("Random seed (0 for unpredictable)")
			;
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

	//measurement method - per timestep or window
	//1D hist or KSG
	//ofile

	fmt::print("Git Head: {}\n", GIT_HEAD);
	fmt::print("Git Log: {}\n", GIT_LOG);
	fmt::print("Git Diff: {}\n", GIT_DIFF);

	fmt::print("Initialising model\n");

	if(target_flock_size > N) PEEXIT("Target flock size is larger than N");
	if(break_apart_threshold >= 1 || break_apart_threshold <= 0) PEEXIT("Break apart threshold should be between 0 and 1");

	std::vector<double> eta_steps;
	{
		std::ifstream eta_fstream(eta_file);
		std::string line;
		while(std::getline(eta_fstream, line))
		{
			eta_steps.push_back(atof(line.c_str()));
		}
		eta_fstream.close();
	}

	seed = calcSeed(seed);
	std::mt19937_64 rng(seed);
	fmt::print("seed : {}\n", seed);

	std::mt19937_64 treeRng(seed);

	FILE* gpipe = popen(gpcmd.c_str(), "w");
	if(gpipe == NULL)                     PEEXIT("failed to open pipe to '{}'\n", gpcmd.c_str());
	if(setvbuf(gpipe, NULL, _IOLBF, 0) != 0) PEEXIT("failed to line-buffer pipe to '{}'\n", gpcmd.c_str());

	const double L = sqrt(((double)N) / rho);
	const double arrowlen = 1;// 0.03*L;

	FlockState main_state(N);
	DoubleBuffer<double> sBuf(N);
	double* const hAgg = (double* const)calloc(N, sizeof(double));

	initialise(imode, N, L, main_state.h.writer(), main_state.x.writer(), 
		main_state.y.writer(), sBuf.writer(), rng);
	main_state.swapBuffers(); sBuf.swap();

	//VP Tree for performing model updates
	VpTreeEuclidean2D modelTree(L, &treeRng);

	// Structure for calculating connected sets

	double rotation = 0;

	//Simulate until one flock reaches target size
	int flockSize = (int)N, flockId = -1;
	size_t targetE = 0;
	bool targetFlockSizeReached = false;
	fmt::print("Finding single flock\n");

	ConnectedSets* connectedSet = new ConnectedSets(N);
	int frameNum = 0;
	

	std::vector<std::vector<double>> hFullData(U_form);
	std::vector<std::vector<double>> hAggFullData(U_form);
	std::vector<std::vector<double>> hNextFullData(U_form);

	for(size_t i = 0; i < U_form; ++i) {
		hFullData[i].resize(N);
		hAggFullData[i].resize(N);
		hNextFullData[i].resize(N);
	}
	std::vector<double> fullGte;
	if(load_state.empty()) {
		int* const rc = (int* const)calloc(N, sizeof(int));
		int* const gc = (int* const)calloc(N, sizeof(int));
		int* const bc = (int* const)calloc(N, sizeof(int));
		for(size_t e = start_eta; e < eta_steps.size(); ++e) {
			fmt::print("Processing eta: {:.2}\n", eta_steps[e]);
			for(size_t u = 0; u < U_form; ++u) {
				connectedSet->clear();
				if(imethod == 0) {
					vicsekMetric(N, L, v, eta_steps[e], main_state, hAgg, rng,
									umethod, modelTree, 0, 0, 0, connectedSet, rotate_frame != 0, &rotation);
				}
				else if(imethod == 1) {
					vicsekTopological(N, L, v, eta_steps[e], dt, main_state, hAgg,
									  rng, topo_neighbours, modelTree, umethod, 0, 0);
				}
				else
					PEEXIT("Invalid imethod");

				auto h = main_state.h.reader();
				auto hnew = main_state.h.writer();
				for(size_t i = 0; i < N; ++i) {
					hFullData[u][i] = h[i];
					hAggFullData[u][i] = hAgg[i];
					hNextFullData[u][i] = hnew[i];
				}

				flockId = getLargestFlock(*connectedSet, flockSize);

				if(render_formation_frames != 0 && u % render_formation_frames == 0) {
					std::string outputDir = "frames/formation";
					make_dir(outputDir);
					fmt::print(gpipe, "set term png size 1280,750\n");
					std::string frameFileName = fmt::format("{}/frame{:04}.png", outputDir, frameNum);
					frameNum++;
					fmt::print(gpipe, "set out \"{}\"\n", frameFileName);
					fmt::print(gpipe, "set termopt enhanced\n");
					//double gteVal = 0;
					fmt::print(gpipe, "set title '{}, N={}, rho={} (L={}), v={}, eta={:.2}, t={}, window={}, largest flock={}, target={}'\n",
						(imethod == 1 ? "Topological" : "Metric"), N, rho, L, v, eta_steps[e], u, windowSize, flockSize, target_flock_size);

					memset(rc, 255, sizeof(int)*N);
					memset(gc, 0, sizeof(int)*N);
					memset(bc, 0, sizeof(int)*N);

					for(size_t i = 0; i < connectedSet->sets[flockId].size(); ++i) {
						rc[connectedSet->sets[flockId][i]] = 0;
						bc[connectedSet->sets[flockId][i]] = 255;
					}
					//fmt::print(gpipe, "set title 'Vicsek Model, N={}'\n", N_single);
					fmt::print(gpipe, "set size square\n");
					fmt::print(gpipe, "set xr [0:{:.16}]\n", L);
					fmt::print(gpipe, "set yr [0:{:.16}]\n", L);
					fmt::print(gpipe, "unset xlabel\n");
					fmt::print(gpipe, "unset ylabel\n");
					fmt::print(gpipe, "set border\n");
					fmt::print(gpipe, "set tics\n");
					fmt::print(gpipe, "rgb(r,g,b) = int(r)*65536 + int(g)*256 + int(b)\n");
					fmt::print(gpipe, "plot '-' u ($2):($3):({0}*cos($1)):({0}*sin($1)):(rgb($4,$5,$6)) w vectors notitle lc rgb variable\n", arrowlen);
					fprintv(gpipe, N, h, main_state.x.reader(), main_state.y.reader(), rc, gc, bc);
					fmt::print(gpipe, "e\n");
					fmt::print(gpipe, "unset out\n");
				}

				if(eta_steps[e] <= eta_threshold && flockSize > (int)target_flock_size) {
					fmt::print("Flock {} has reached size of {} members at t={} of \\eta={:.2}\n", flockId, flockSize, u, eta_steps[e]);
					targetE = e;
					targetFlockSizeReached = true;
					break;
				}

				// swap buffers
				main_state.swapBuffers(); sBuf.swap();

				progrep("sim ", u, U_form);

			}

			if(calc_cooling_gte == 1) {
				DataSet<3> gteData((int)(U_form * N), (int)KSG_neighbours);

				int j = 0;
				for(size_t u1 = 0; u1 < U_form; ++u1) {
					for(size_t i = 0; i < N; ++i, ++j) {
						gteData.setData(hFullData[u1][i], j, 0);
						gteData.setData(hAggFullData[u1][i], j, 1);
						gteData.setData(hNextFullData[u1][i], j, 2);
					}
				}

				fullGte.push_back(calculateTE(gteData));
			}

			if(targetFlockSizeReached)
				break;
		}
		if(calc_cooling_gte == 1) {
			FILE * gteFP = fopen("fullgte.bin", "wb");
			size_t numgte = fullGte.size();
			if(fwrite(&numgte, sizeof(size_t), 1, gteFP) != 1) PEEXIT("write failed");
			for(std::vector<double>::reverse_iterator itr = fullGte.rbegin(); itr != fullGte.rend(); itr++) {
				double gteVal = *itr;
				if(fwrite(&gteVal, sizeof(double), 1, gteFP) != 1) PEEXIT("write failed");
			}
			fclose(gteFP);
		}

		free(rc); free(gc); free(bc);
	}
	else {
		fmt::print("Loading state from file {}\n", load_state);
		FILE* loadFP = fopen(load_state.c_str(), "rb");
		if(fread(main_state.h.writer(), sizeof(double), N, loadFP) != N) PEEXIT("load failed");
		if(fread(main_state.x.writer(), sizeof(double), N, loadFP) != N) PEEXIT("load failed");
		if(fread(main_state.y.writer(), sizeof(double), N, loadFP) != N) PEEXIT("load failed");
		main_state.swapBuffers();
		if(fread(main_state.h.writer(), sizeof(double), N, loadFP) != N) PEEXIT("load failed");
		if(fread(main_state.x.writer(), sizeof(double), N, loadFP) != N) PEEXIT("load failed");
		if(fread(main_state.y.writer(), sizeof(double), N, loadFP) != N) PEEXIT("load failed");
		if(fread(hAgg, sizeof(double), N, loadFP) != N) PEEXIT("load failed");
		load(*connectedSet, loadFP);
		fclose(loadFP);
		flockId = getLargestFlock(*connectedSet, flockSize);

	}

	if(!save_state.empty()) {
		fmt::print("Saving state to file {}\n", save_state);
		FILE* saveFP = fopen(save_state.c_str(), "wb");
		if(fwrite(main_state.h.reader(), sizeof(double), N, saveFP) != N) PEEXIT("save failed");
		if(fwrite(main_state.x.reader(), sizeof(double), N, saveFP) != N) PEEXIT("save failed");
		if(fwrite(main_state.y.reader(), sizeof(double), N, saveFP) != N) PEEXIT("save failed");
		if(fwrite(main_state.h.writer(), sizeof(double), N, saveFP) != N) PEEXIT("save failed");
		if(fwrite(main_state.x.writer(), sizeof(double), N, saveFP) != N) PEEXIT("save failed");
		if(fwrite(main_state.y.writer(), sizeof(double), N, saveFP) != N) PEEXIT("save failed");
		if(fwrite(hAgg, sizeof(double), N, saveFP) != N) PEEXIT("save failed");
		serialize(*connectedSet, saveFP);
		fclose(saveFP);
	}


	if(!targetFlockSizeReached) {
		fmt::print("Target flock size ({}) was not reached after last \\eta step. Continuing with largest flock at last time step. Flock size is {}\n", target_flock_size, flockSize);
		targetE = eta_steps.size() - 1;
	}

	const size_t N_single = dont_cull ? N : flockSize;// flockSize;
	FlockState saveState(N_single);
	double* const hAgg_saveState = (double* const)calloc(N_single, sizeof(double));

	//Save the state with the target flock
	//Use a temporary FlockState to have two writers, feed "h,x,y" values into writer of
	//first, and "hnew,xnew,ynew" into writer of second. Swap first, and then copy from
	//writer of second into writer of first.
	fmt::print("Saving flock state\n");
	if(!dont_cull && flockId != -1) {
		FlockState saveState_b(N_single);
		auto h = main_state.h.reader(), x = main_state.x.reader(), y = main_state.y.reader();
		auto hnew = main_state.h.writer(), xnew = main_state.y.writer(), ynew = main_state.y.writer();

		auto h_ss = saveState.h.writer(), x_ss = saveState.x.writer(), y_ss = saveState.y.writer();
		auto hnew_ss = saveState_b.h.writer(), xnew_ss = saveState_b.x.writer(), ynew_ss = saveState_b.y.writer();
		for(size_t i = 0, j = 0; i < connectedSet->sets[flockId].size(); ++i) {
			int pId = connectedSet->sets[flockId][i];
			h_ss[j] = h[pId];
			x_ss[j] = x[pId];
			y_ss[j] = y[pId];
			hnew_ss[j] = hnew[pId];
			xnew_ss[j] = xnew[pId];
			ynew_ss[j] = ynew[pId];
			hAgg_saveState[j] = hAgg[pId];
			++j;
		}
		saveState.swapBuffers();
		std::copy(hnew_ss, hnew_ss + N_single, saveState.h.writer());
		std::copy(xnew_ss, xnew_ss + N_single, saveState.x.writer());
		std::copy(ynew_ss, ynew_ss + N_single, saveState.y.writer());
	}
	else {
		FlockState saveState_b(N_single);
		auto h = main_state.h.reader(), x = main_state.x.reader(), y = main_state.y.reader();
		auto hnew = main_state.h.writer(), xnew = main_state.y.writer(), ynew = main_state.y.writer();

		auto h_ss = saveState.h.writer(), x_ss = saveState.x.writer(), y_ss = saveState.y.writer();
		auto hnew_ss = saveState_b.h.writer(), xnew_ss = saveState_b.x.writer(), ynew_ss = saveState_b.y.writer();

		for(size_t i = 0, j = 0; i < N_single; ++i) {
			int pId = (int)i;
			h_ss[j] = h[pId];
			x_ss[j] = x[pId];
			y_ss[j] = y[pId];
			hnew_ss[j] = hnew[pId];
			xnew_ss[j] = xnew[pId];
			ynew_ss[j] = ynew[pId];
			hAgg_saveState[j] = hAgg[pId];
			++j;
		}
		saveState.swapBuffers();
		std::copy(hnew_ss, hnew_ss + N_single, saveState.h.writer());
		std::copy(xnew_ss, xnew_ss + N_single, saveState.x.writer());
		std::copy(ynew_ss, ynew_ss + N_single, saveState.y.writer());
	}

	delete connectedSet;

	free(hAgg);

	//Adjust buffers to new size
	FlockState fs_single(N_single);
	double* const hAgg_single = (double* const)calloc(N_single, sizeof(double));

	VpTreeEuclidean2D newModelTree(L, &treeRng);
	ConnectedSets smallSet(N_single);
	fmt::print("Breaking flock apart\n");
	const int numEta = (int)targetE + 1;
	std::vector<std::vector<double>> gte(numEta);
	std::vector<std::vector<double>> miFlock(numEta);
	std::vector<std::vector<double>> alpha(numEta);
	std::vector<std::vector<double>> order(numEta);
	std::vector<std::vector<size_t>> flockSizes(numEta);
	std::vector<size_t> breakApartTime(numEta);
	for(int i = 0; i < numEta; ++i) {
		breakApartTime[i] = 0;
	}

	std::vector<std::vector<double>> hData(windowSize); //h_t for particles in flock at time t
	std::vector<std::vector<double>> hNextData(windowSize); //h_{t+1} for particles in flock at time t
	std::vector<std::vector<double>> hAggData(windowSize); //hAgg_t for particles in flock at time t
	std::vector<double> hAvgData(windowSize);

	std::vector<int> rVec, gVec, bVec;

	//Simulate the single flock at all higher noise values for time or until flock breaks up. Track GTE
	for(int e = (int)targetE; e >= 0; --e) {
		//Load the saved state
		fs_single = saveState;
		std::copy(hAgg_saveState, hAgg_saveState + N_single, hAgg_single);
		
		fmt::print("Trying \\eta[{}]={:.2}\n", e, eta_steps[e]);
		int frameCount = 0;
		for(size_t u = 0; u < U_break; ++u) {
			auto h = fs_single.h.reader(), x = fs_single.x.reader(), y = fs_single.y.reader();
			auto hnew = fs_single.h.writer();
			smallSet.clear();
			

			if(imethod == 0) {
				vicsekMetric(N_single, L, v, eta_steps[e], fs_single, hAgg_single, rng,
								umethod, newModelTree, 0, 0, 0, &smallSet, rotate_frame != 0, &rotation);
			}
			else if(imethod == 1) {
				vicsekTopological(N_single, L, v, eta_steps[e], dt, fs_single, hAgg_single,
								  rng, topo_neighbours, newModelTree, umethod, 0, 0);
			}
			else
				PEEXIT("Invalid imethod");
			double zx = 0.0, zy = 0.0;
			if(!calc_order_just_lcc) {
				int count = 0;
				for(size_t i = 0; i < smallSet.particleMap.size(); ++i) {
					if(smallSet.particleMap[i] == flockId) {
						double ux, uy;
						sincos(h[i], &uy, &ux); // u is unit vector in direction h[i]
						zx += ux;
						zy += uy;
						++count;
					}
				}
				zx /= (double)count;
				zy /= (double)count;
			}
			else {
				order_param(N_single, h, &zx, &zy);
			}
			order[e].push_back(sqrt(zx*zx + zy*zy));

			
			flockId = getLargestFlock(smallSet, flockSize);

			double xavg = 0.0, yavg = 0.0;
			for(size_t i = 0; i < smallSet.particleMap.size(); ++i) {
				if(smallSet.particleMap[i] == flockId) {
					double xangle, yangle;
					sincos(h[i], &yangle, &xangle);
					xavg += xangle;
					yavg += yangle;
				}
			}

			//Copy all the data. Terry also mentioned using just the data for the particles in the flock
			size_t curIdx = u % windowSize;
			hData[curIdx].clear();
			hNextData[curIdx].clear();
			hAggData[curIdx].clear();
			hAvgData[curIdx] = awrap(atan2(yavg, xavg));
			rVec.clear();
			gVec.clear();
			bVec.clear();
			if(calc_metric_allparticles == 1) {
				for(size_t i = 0; i < N_single; ++i) {
					hData[curIdx].push_back(h[i]);
					hNextData[curIdx].push_back(hnew[i]);
					hAggData[curIdx].push_back(hAgg[i]);
					rVec.push_back(0); gVec.push_back(0); bVec.push_back(255);
				}
			}
			else {
				for(size_t i = 0; i < smallSet.particleMap.size(); ++i) {
					if(smallSet.particleMap[i] == flockId) {
						hData[curIdx].push_back(h[i]);
						hNextData[curIdx].push_back(hnew[i]);
						hAggData[curIdx].push_back(hAgg[i]);
						rVec.push_back(0); gVec.push_back(0); bVec.push_back(255);
					}
					else {
						rVec.push_back(255); gVec.push_back(0); bVec.push_back(0);
					}
				}
			}

			if(calc_mi_against_flock == 1 && u > 5) {
				int numFrames = std::min((int)u, windowSize);
				int numElements = 0;
				for(int i = 0; i < numFrames; ++i) {
					numElements += (int)hData[i].size();
				}
				if(numElements > 0) {
					DataSet<2> miData(numElements, (int)KSG_neighbours);

					int j = 0;
					for(size_t u1 = 0; u1 < (size_t)numFrames; ++u1) {
						for(size_t i = 0; i < hData[u1].size(); ++i, ++j) {
							miData.setData(hData[u1][i], j, 0);
							miData.setData(hAvgData[u1], j, 1);
						}
					}

					miFlock[e].push_back(calculateMI(miData));
					if(!miFlock[e].empty() && !order[e].empty()) {
						alpha[e].push_back(miFlock[e].back() / order[e].back());
					}
				}
			}

			//Calc GTE over the last windowSize timesteps of the flock
			if(calc_gte == 1 && u >= 1) {
				int numFrames = std::min((int)u, windowSize);
				int numElements = 0;// numFrames * N_single;
				for(int i = 0; i < numFrames; ++i) {
					numElements += (int)hData[i].size();
				}
				DataSet<3> gteData(numElements, (int)KSG_neighbours);

				int j = 0;
				for(int u1 = 0; u1 < numFrames; ++u1) {
					for(size_t i = 0; i < hData[u1].size(); ++i, ++j) {
						gteData.setData(hData[u1][i], j, 0);
						gteData.setData(hAggData[u1][i], j, 1);
						gteData.setData(hNextData[u1][i], j, 2);
					}
				}
				gte[e].push_back(calculateTE(gteData));
				flockSizes[e].push_back(flockSize);
			}

			if(render_break_frames != 0 && u % render_break_frames == 0) {
				std::string outputDir = fmt::format("frames/eta_{:.2}", eta_steps[e]);
				make_dir(outputDir);
				//fmt::print(gpipe, "set terminal postscript eps enhanced size 6in,6in\n");
				fmt::print(gpipe, "set term png size 1280,750\n");
				std::string frameFileName = fmt::format("{}/frame{:04}.png", outputDir, frameCount);
				frameCount++;
				fmt::print(gpipe, "set out \"{}\"\n", frameFileName);
				fmt::print(gpipe, "set termopt enhanced\n");
				double gteVal = 0;
				if(!gte[e].empty())
					gteVal = gte[e].back();
				fmt::print(gpipe, "set multiplot layout 1,3 title '{}, N={}({}), rho={} (L={}), v={}, eta={:.2}f, GTE={:.2} Nats, t={}, window={}, largest flock={}, target={}'\n",
					(imethod == 1 ? "Topological" : "Metric"), N, N_single, rho, L, v, eta_steps[e], gteVal, u, windowSize, flockSize, target_flock_size);

				fmt::print(gpipe, "set title 'Flock'\n");
				fmt::print(gpipe, "set size square\n");
				fmt::print(gpipe, "set xr [0:{:.16}]\n", L);
				fmt::print(gpipe, "set yr [0:{:.16}]\n", L);
				fmt::print(gpipe, "unset xlabel\n");
				fmt::print(gpipe, "unset ylabel\n");
				fmt::print(gpipe, "set border\n");
				fmt::print(gpipe, "set tics\n");
				fmt::print(gpipe, "rgb(r,g,b) = int(r)*65536 + int(g)*256 + int(b)\n");
				fmt::print(gpipe, "plot '-' u ($2):($3):({0}*cos($1)):({0}*sin($1)):(rgb($4,$5,$6)) w vectors notitle lc rgb variable\n", arrowlen);
				fprintv(gpipe, N_single, h, x, y, &(rVec[0]), &(gVec[0]), &(bVec[0]));
				fmt::print(gpipe, "e\n");

				fmt::print(gpipe, "set title 'GTE of blue particles'\n");
				fmt::print(gpipe, "set xr [1:{}]\n", ((u / 100)+1) * 100);
				fmt::print(gpipe, "set yr [0:1.1]\n");
				fmt::print(gpipe, "set xlabel 'Time'\n");
				fmt::print(gpipe, "unset ylabel\n");
				fmt::print(gpipe, "set border\n");
				fmt::print(gpipe, "set tics\n");
				fmt::print(gpipe, "set key\n");
				fmt::print(gpipe, "plot");
				int plotCount = 0;
				if(!gte[e].empty()) {
					fmt::print(gpipe, " '-' u 1:2 w lines title 'GTE' lc rgb \"red\"");
					++plotCount;
				}

				if(!order[e].empty()) {
					if(plotCount != 0)
						fmt::print(gpipe, ",");
					fmt::print(gpipe, " '-' u 1:2 w lines title 'Order' lc rgb \"blue\"");
					++plotCount;
				}

				if(!miFlock[e].empty()) {
					if(plotCount != 0)
						fmt::print(gpipe, ",");
					fmt::print(gpipe, " '-' u 1:2 w lines title 'MI against Flock' lc rgb \"green\"");
					++plotCount;
				}

				if(plotCount == 0) {
					fmt::print(gpipe, " 2 with lines");
				}
				fmt::print(gpipe, "\n");

				if(!gte[e].empty()) {
					for(size_t i = 0; i < gte[e].size(); ++i) {
						fmt::print(gpipe, "{:20.16} {:20.16}\n", (float)i, gte[e][i]);
					}
					fmt::print(gpipe, "e\n");
				}
				if(!order[e].empty()) {
					for(size_t i = 0; i < order[e].size(); ++i) {
						fmt::print(gpipe, "{:20.16} {:20.16}\n", (float)i, order[e][i]);
					}
					fmt::print(gpipe, "e\n");
				}
				if(!miFlock[e].empty()) {
					for(size_t i = 0; i < miFlock[e].size(); ++i) {
						fmt::print(gpipe, "{:20.16} {:20.16}\n", (float)i, miFlock[e][i]);
					}
					fmt::print(gpipe, "e\n");
				}
				if(plotCount == 0) {
					fmt::print(gpipe, "e\n");
				}

				double minAlpha = 0., maxAlpha = 1.;
				if(!alpha[e].empty()) {
					minAlpha = *std::min_element(alpha[e].begin(), alpha[e].end());
					maxAlpha = *std::max_element(alpha[e].begin(), alpha[e].end());
				}
				double avgAlpha = 0.0;
				for(size_t i = 0; i < alpha[e].size(); ++i) {
					avgAlpha += alpha[e][i];
				}
				avgAlpha /= (double)alpha[e].size();
				maxAlpha += (minAlpha == maxAlpha); //If equal, add one to maxAlpha
				fmt::print(gpipe, "set title 'MI / Order - Min: {}, Max: {}, Avg: {}'\n", minAlpha, maxAlpha, avgAlpha);
				fmt::print(gpipe, "set xr [1:{}]\n", ((u / 100) + 1) * 100);
				fmt::print(gpipe, "set yr [{}:{}]\n", minAlpha, maxAlpha);
				fmt::print(gpipe, "set xlabel 'Time'\n");
				fmt::print(gpipe, "unset ylabel\n");
				fmt::print(gpipe, "set border\n");
				fmt::print(gpipe, "set tics\n");
				fmt::print(gpipe, "set key\n");
				fmt::print(gpipe, "plot '-' u 1:2 w lines title 'Alpha' lc rgb \"red\"\n");
				if(!alpha[e].empty()) {
					for(size_t i = 0; i < alpha[e].size(); ++i) {
						fmt::print(gpipe, "{:20.16} {:20.16}\n", (float)i, alpha[e][i]);
					}
				}
				else
					fmt::print(gpipe, "0 0\n");
				fmt::print(gpipe, "e\n");

				fmt::print(gpipe, "unset multiplot\n");
				fmt::print(gpipe, "unset out\n");
			}


			if(breakApartTime[e] == 0 && flockSize < (int)((double)target_flock_size * break_apart_threshold)) {
				fmt::print("Flock {} has broken apart ({} members left) at t={} of \\eta={:.2}\n", flockId, flockSize, u, eta_steps[e]);
				breakApartTime[e] = u;
				//break;
			}

			// swap buffers
			fs_single.swapBuffers();

			progrep("sim ", u, U_break);
		}
	}
	pclose(gpipe);

	//Write GTE info to file
	FILE* fp = fopen("singleflock.bin", "wb");
	if(fwrite(&numEta, sizeof(int), 1, fp) != 1) PEEXIT("write failed");
	for(int i = 0; i < numEta; ++i) {
		double etaVal = eta_steps[i];
		size_t breakTime = breakApartTime[i];
		int etaCount = (int)gte[i].size();
		if(fwrite(&etaVal, sizeof(double), 1, fp) != 1) PEEXIT("write failed");
		if(fwrite(&breakTime, sizeof(size_t), 1, fp) != 1) PEEXIT("write failed");
		if(fwrite(&etaCount, sizeof(int), 1, fp) != 1) PEEXIT("write failed");
		if((int)fwrite(&(gte[i][0]), sizeof(double), etaCount, fp) != etaCount) PEEXIT("write failed");
		if((int)fwrite(&(flockSizes[i][0]), sizeof(size_t), etaCount, fp) != etaCount) PEEXIT("write failed");
	}
	fclose(fp);

	free(hAgg_saveState);

	return EXIT_SUCCESS;
}

#include "../include/ksgControl.hpp"
#include "../include/dataSet.hpp"
#include "../include/metrics.hpp"
#include "../include/defs.hpp"
#include "../include/logFileHelper.hpp"
#include "../include/nearestNeighbour.hpp"

#include <fmt/format.h>
#include <stdio.h>

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
							  std::mt19937_64& rng
)
{
	const size_t recordedU = sp.U / sp.record_T_steps;
	auto const h = fh.read_h();
	auto const x = fh.read_x();
	auto const y = fh.read_y();
	ITStats itStats;

	if(sp.sims[MIBIN])
	{
		ScopeTimer st("MI (flock)");

		fmt::print("Calculating MI\n");
		size_t numElements = gp.N * recordedU;
		fmt::print("Num particles: {}\n", (int)numElements);
		
		DataSet<2> dataForMI((int)numElements, (int)sp.KSG_neighbours);
		int j = 0;
		for(size_t u = 0; u < recordedU; ++u) {
			//size_t elements = interactions[u].size() / 2;
			for(size_t i = 0; i < gp.N; i++) {
				//int p1 = interactions[u][i * 2];
				double sub = 0;// avgDirections[u][flockIds[u*p.N + i]];
				dataForMI.setData(awrap(h[u + recordedU * i] - cumulativeRotations[u] - sub), j, 0);
				dataForMI.setData(avgDirections[u][flockIds[u*gp.N + i]], j, 1);

				//int p2 = interactions[u][i * 2 + 1];
				//dataForMI.setData(awrap(hFull[u + p.U * p1] - cumulativeRotations[u] - avgDirections[u][flockIds[u*p.N + p1]]), j, 0);
				//dataForMI.setData(awrap(hFull[u + p.U * p2] - cumulativeRotations[u] - avgDirections[u][flockIds[u*p.N + p2]]), j, 1);
				j++;
			}
		}

		if(sp.shuffle)
		{
			std::vector<double> shuffleVec(numElements);
			size_t shuffle_idx = 0;
			for(size_t u = 0; u < recordedU; ++u) {
				size_t elements = interactions[u].size() / 2;
				for(size_t i = 0; i < elements; i++) {
					shuffleVec[shuffle_idx] = h[u + recordedU * interactions[u][i * 2 + 1]];
					shuffle_idx++;
				}
			}
			std::shuffle(shuffleVec.begin(), shuffleVec.end(), rng);
			dataForMI.setData(&(shuffleVec[0]), 1);
		}
		fmt::print("numElements: {}, j: {}\n", (int)numElements, j);

		/*bool verified = verifyMI(dataForMI, 0.0);
		if(!verified)
		fmt::print("MI not working right\n");*/
		double val = calculateMI(dataForMI);
		itStats.mi.accumulate(val);
	}

	if(sp.sims[TEBIN])
	{
		ScopeTimer st("TE (flock)");

		fmt::print("Calculating TE\n");
		size_t U1 = recordedU - 1; //Reduced by one for tebin and gtebin, for past of X. 

		size_t numTEElements = 0;
		for(size_t i = 0; i < U1; ++i) {
			numTEElements += interactions[i].size() / 2;
		}
		fmt::print("All interactions: {}\n", (int)numTEElements);

		DataSet<3> dataForTE((int)numTEElements, (int)sp.KSG_neighbours);
		int j = 0;
		for(size_t u = 0; u < U1; ++u) {
			size_t numElements = (int)(interactions[u].size()) / 2;
			for(size_t i = 0; i < numElements; i++, j++/*, j += 2*/) {
				int p1 = interactions[u][i * 2];		//First particle in interaction pair
				int p2 = interactions[u][i * 2 + 1];	//Second particle in interaction pair
				int t = (int)(u);		//Current timestep (i.e. for Xt and Yt)
				int t1 = (int)(u + 1);			//Next timestep (i.e. for Xt+1)

				dataForTE.setData(awrap(h[t + p1 *recordedU] - cumulativeRotations[t]), j, 0);
				dataForTE.setData(awrap(h[t + p2 *recordedU] - cumulativeRotations[t]), j, 1);
				dataForTE.setData(awrap(h[t1 + p1 * recordedU] - cumulativeRotations[t]), j, 2);
			}
		}
		if(sp.shuffle)
		{
			std::vector<double> shuffleVec(numTEElements);
			size_t shuffle_idx = 0;
			for(size_t u = 0; u < U1; ++u) {
				int t1 = (int)u;
				if(sp.te_shuffle_dim == 2)
					t1++;
				size_t elements = interactions[u].size() / 2;
				for(size_t i = 0; i < elements; i++) {
					int pId = (int)i * 2;
					if(sp.te_shuffle_dim == 1)
						pId++;
					shuffleVec[shuffle_idx] = h[t1 + recordedU * interactions[u][pId]];
					shuffle_idx++;
				}
			}
			std::shuffle(shuffleVec.begin(), shuffleVec.end(), rng);
			dataForTE.setData(&(shuffleVec[0]), sp.te_shuffle_dim);
		}
		fmt::print("numElements: {}, j: {}\n", (int)numTEElements, j);
		//bool res = verifyTE(dataForTE, 0.0);
		//fmt::print("Verify result: {}\n", res ? "success" : "failure");
		double val = calculateTE(dataForTE);
		itStats.te.accumulate(val);
	}

	if(sp.sims[GTEBIN])
	{
		ScopeTimer st("GTE (flock)");

		fmt::print("Calculating GTE\n");
		size_t U1 = recordedU - 1; //Reduced by one for tebin and gtebin, for past of X. 
		double gte = 0.0;

		if(sp.ksg_gte_dims == 3) {
			DataSet<3> dataForGTE((int)(gp.N * (recordedU - 1)), (int)sp.KSG_neighbours);

			int j = 0;
			for(size_t u = 0; u < U1; ++u) {
				for(size_t i = 0; i < gp.N; i++, ++j) {
					int t = (int)(u);		//Current timestep (i.e. for Xt and Yt)
					int t1 = (int)(u + 1);			//Next timestep (i.e. for Xt+1)

					dataForGTE.setData(awrap(h[t + i * recordedU] - cumulativeRotations[t]), j, 0);
					dataForGTE.setData(awrap(hAggFull[t + i * recordedU] - cumulativeRotations[t]), j, 1);
					dataForGTE.setData(awrap(h[t1 + i * recordedU] - cumulativeRotations[t]), j, 2);
				}
			}

			if(sp.shuffle)
			{
				std::vector<double> shuffleVec(gp.N * (recordedU - 1));
				size_t shuffle_idx = 0;
				for(size_t u = 0; u < U1; ++u) {
					int t1 = (int)u;
					if(sp.te_shuffle_dim == 2)
						t1++;
					for(size_t i = 0; i < gp.N; i++) {
						if(sp.te_shuffle_dim == 1)
							shuffleVec[shuffle_idx] = hAggFull[t1 + recordedU * i];
						else
							shuffleVec[shuffle_idx] = h[t1 + recordedU * i];
						shuffle_idx++;
					}
				}
				std::shuffle(shuffleVec.begin(), shuffleVec.end(), rng);
				dataForGTE.setData(&(shuffleVec[0]), sp.te_shuffle_dim);
			}
			fmt::print("numElements: {}, j: {}\n", (int)(gp.N*(recordedU - 1)), j);
			gte = calculateTE(dataForGTE);
		}
		else if(sp.ksg_gte_dims == 0) {
			DataSet<GK + 2> dataForGTE((int)(gp.N * (recordedU - 1)), (int)sp.KSG_neighbours);

			int j = 0;
			for(size_t u = 0; u < U1; ++u) {
				const size_t neigh_t = u; //Timestep to check for topo_neighbours nearest neighbours
				for(size_t l = 0; l < gp.N; ++l) {
					posPts[l].values[0] = x[neigh_t + l*recordedU];
					posPts[l].values[1] = y[neigh_t + l*recordedU];
					posPts[l].idx = (int)l;
				}
				getNearestNeighbours(posPts, (int)gp.N, (int)GK, gp.L, fullGTENNidx);

				std::vector<double> neighbourAngles(GK);
				for(size_t i = 0; i < gp.N; i++, ++j) {
					int t = (int)(u);		//Current timestep (i.e. for Xt and Yt)
					int t1 = (int)(u + 1);			//Next timestep (i.e. for Xt+1)

					dataForGTE.setData(h[t + i * recordedU], j, 0);
					for(size_t y_i = 0; y_i < (size_t)GK; ++y_i) {
						neighbourAngles[y_i] = h[t + fullGTENNidx[i*GK + y_i] * recordedU];
					}
					//std::sort(neighbourAngles.begin(), neighbourAngles.end());
					for(size_t y_i = 0; y_i < (size_t)GK; ++y_i) {
						dataForGTE.setData(neighbourAngles[y_i], j, (int)y_i + 1);
					}
					dataForGTE.setData(h[t1 + i * recordedU], j, GK + 1);
				}
			}
			fmt::print("numElements: {}, j: {}\n", (int)(gp.N*(recordedU - 1)), j);
			gte = calculateGTE(dataForGTE);
		}
		else {
			PEEXIT("Flock Particle KSG not supported for ksg_gte_dims={}", sp.ksg_gte_dims);
		}
		itStats.gte.accumulate(gte);
	}
	return itStats;
}

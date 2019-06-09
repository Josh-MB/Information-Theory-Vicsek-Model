#include "../include/ksgControl.hpp"
#include "../include/dataSet.hpp"
#include "../include/metrics.hpp"
#include "../include/defs.hpp"
#include "../include/logFileHelper.hpp"
#include "../include/nearestNeighbour.hpp"

#include <fmt/format.h>
#include <stdio.h>

ITStats fillAndRunSubsetKSG(GeneratorParams const& gp,
						 SimParams const& sp,
						 std::vector<std::vector<int> > const& interactions,
						 FlockHistory const& fh,
						 double const*const hAggFull,
						 double const*const cumulativeRotations,
						 std::mt19937_64& rng
						)
{
	fmt::print("Using subsets. Subset_size: {}\n", sp.subset_size);
	const size_t recordedU = sp.U / sp.record_T_steps;
	int subset_runs_use = std::max(1, sp.subset_runs);
	auto const h = fh.read_h();
	//if(subset_size_use > numElements * 0.9)
	//	PEEXIT("Subset size was more than 90% of points: {} ({})", subset_size_use, numElements);
	ITStats itStats;

	fmt::print("Calculating MI and TE (if requested) using subsets\n");
	{
		size_t numElements = 0;
		for(size_t i = 0; i < recordedU; ++i) {
			numElements += interactions[i].size();
		}
		numElements /= 2;

		int subset_size_use = 1000;
		if(sp.subset_size <= 1.0)
			subset_size_use = (int)((double)numElements * sp.subset_size);
		else
			subset_size_use = std::max(1000, (int)sp.subset_size);

		fmt::print("NumElements: {}, Using {} subsets with window size: {}\n", numElements, subset_runs_use, subset_size_use);

		std::vector<InteractionPair> interactionIndices;
		interactionIndices.reserve(numElements);

		for(size_t i = 0; i < recordedU; ++i) {
			int num = (int)interactions[i].size() / 2;
			for(int j = 0; j < num; ++j) {
				interactionIndices.push_back(InteractionPair((int)i, j));
			}
		}
		fmt::print("Interaction indices filled\n");

		for(int run = 0; run < subset_runs_use; ++run) {
			//For MI and TE
			std::vector<ksg_point> subset_points(subset_size_use);

			fmt::print("Progress for run {}:\n", run);

			shuffleFirstMIndices(interactionIndices, subset_size_use, rng);

			for(int it = 0; it < subset_size_use; it++) {
				int randu = interactionIndices[it].first;
				int randi = interactionIndices[it].second;

				subset_points[it].h1 = awrap(h[randu + recordedU * interactions[randu][randi * 2]] - cumulativeRotations[randu]);
				subset_points[it].h2 = awrap(h[randu + recordedU * interactions[randu][randi * 2 + 1]] - cumulativeRotations[randu]);
				subset_points[it].h1new = awrap(h[randu + 1 + recordedU * interactions[randu][randi * 2]] - cumulativeRotations[randu]);
				subset_points[it].h2Agg = awrap(hAggFull[randu + recordedU * interactions[randu][randi * 2]] - cumulativeRotations[randu]);
			}

			fmt::print("\tPoints chosen for MI and TE...\n");

			if(sp.sims[MIBIN])
			{
				DataSet<2> dataForMI(subset_size_use, (int)sp.KSG_neighbours);
				for(int i = 0; i < subset_size_use; ++i)
				{
					dataForMI.setData(subset_points[i].h1, i, 0);
					dataForMI.setData(subset_points[i].h2, i, 1);
				}

				double val = calculateMI(dataForMI);
				itStats.mi.accumulate(val);
			}

			fmt::print("\tMI calculated...\n");

			if(sp.sims[TEBIN])
			{
				DataSet<3> dataForTE(subset_size_use, (int)sp.KSG_neighbours);
				for(int i = 0; i < subset_size_use; ++i)
				{
					dataForTE.setData(subset_points[i].h1, i, 0);
					dataForTE.setData(subset_points[i].h2, i, 1);
					dataForTE.setData(subset_points[i].h1new, i, 2);
				}

				double val = calculateTE(dataForTE);
				itStats.te.accumulate(val);
			}
			fmt::print("\tTE calculated...\n");
			fmt::print("\tSubset {} done\n", run);
		}
	}

	fmt::print("Calculating GTE (if requested) using subsets\n");
	{
		int subset_size_use = 1000;
		if(sp.subset_size <= 1.0)
			subset_size_use = (int)((double)(gp.N*recordedU) * sp.subset_size);
		else
			subset_size_use = std::max(1000, (int)sp.subset_size);

		fmt::print("NumElements: {}, Using {} subsets with window size: {}\n", gp.N*recordedU, subset_runs_use, subset_size_use);

		std::vector<InteractionPair> indices;
		indices.reserve(gp.N*recordedU);

		for(size_t i = 0; i < recordedU; ++i) {
			for(int j = 0; j < (int)gp.N; ++j) {
				indices.push_back(InteractionPair((int)i, j));
			}
		}
		fmt::print("Indices filled\n");

		for(int run = 0; run < subset_runs_use; ++run) {

			std::vector<ksg_point> subset_points(subset_size_use);

			fmt::print("Progress for run {}:\n", run);

			shuffleFirstMIndices(indices, subset_size_use, rng);

			for(int it = 0; it < subset_size_use; it++) {
				int randu = indices[it].first;
				int randi = indices[it].second;

				subset_points[it].h1 = awrap(h[randu + recordedU * randi] - cumulativeRotations[randu]);
				subset_points[it].h1new = awrap(h[randu + 1 + recordedU * randi] - cumulativeRotations[randu]);
				subset_points[it].h2Agg = awrap(hAggFull[randu + recordedU * randi] - cumulativeRotations[randu]);
			}

			if(sp.sims[GTEBIN])
			{
				DataSet<3> dataForGTE(subset_size_use, (int)sp.KSG_neighbours);
				for(int i = 0; i < subset_size_use; ++i)
				{
					dataForGTE.setData(subset_points[i].h1, i, 0);
					dataForGTE.setData(subset_points[i].h2Agg, i, 1);
					dataForGTE.setData(subset_points[i].h1new, i, 2);
				}

				double val = calculateTE(dataForGTE);
				itStats.gte.accumulate(val);
			}
			fmt::print("\tGTE calculated...\n");
			fmt::print("\tSubset {} done\n", run);
		}
	}
	return itStats;
}

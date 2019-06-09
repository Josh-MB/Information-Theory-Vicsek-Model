#include "../include/ksgControl.hpp"
#include "../include/dataSet.hpp"
#include "../include/metrics.hpp"
#include "../include/defs.hpp"
#include "../include/nearestNeighbour.hpp"

#include <fmt/format.h>
#include <stdio.h>

double fillAndRunPairwiseMI(std::vector<int> const & interactions, double const * const h, size_t const KSG_neighbours)
{
	int numElements = (int)interactions.size() / 2;
	DataSet<2> dataForMI((int)numElements, (int)KSG_neighbours);
	for(int i = 0; i < numElements; i++) {
		int p1 = interactions[i * 2];
		int p2 = interactions[i * 2 + 1];
		//dataForMI.setData(count[u + U * p1], i, 0);
		//dataForMI.setData(count[u + U * p2], i, 1);
		dataForMI.setData(h[p1], i, 0);
		dataForMI.setData(h[p2], i, 1);
	}

	return calculateMI(dataForMI);
}

double fillAndRunPairwiseTE(std::vector<int> const & interactions, double const * const h, double const * const hnew, size_t const KSG_neighbours)
{
	const int numElements = (int)interactions.size() / 2;
	//const size_t numTEElements = numElements *2;

	DataSet<3> dataForTE((int)numElements, (int)KSG_neighbours);

	int j = 0;
	for(int i = 0; i < numElements; i++, j++/*, j+=2*/) {
		int p1 = interactions[i * 2];		//First particle in interaction pair
		int p2 = interactions[i * 2 + 1];	//Second particle in interaction pair

		dataForTE.setData(h[p1], j, 0);
		dataForTE.setData(h[p2], j, 1);
		dataForTE.setData(hnew[p1], j, 2);
		/*dataForTE.setData(h[p2], j + 1, 0);
		dataForTE.setData(h[p1], j + 1, 1);
		dataForTE.setData(hnew[p2], j + 1, 2);*/
	}
	fmt::print("Dataset points placed: {}, numElements: {}\n", (int)j, (int)numElements);

	return calculateTE(dataForTE);
}

double fillAndRunPairwiseGTE(GeneratorParams const& gp,
							 SimParams const& sp,
							 std::vector<VPoint<2> > & posPts,
							 double const*const x,
							 double const*const y,
							 double const*const h,
							 double const*const hnew,
							 double const*const hAgg,
							 int * fullGTENNidx)
{
	if(sp.ksg_gte_dims == 3) {
		DataSet<3> dataForGTE((int)gp.N, (int)sp.KSG_neighbours);

		for(int i = 0; i < (int)gp.N; ++i) {
			dataForGTE.setData(h[i], i, 0);
			dataForGTE.setData(hAgg[i], i, 1);
			dataForGTE.setData(hnew[i], i, 2);
		}
		return calculateTE(dataForGTE);
	}
	else if(sp.ksg_gte_dims==0) {
		for(int j = 0; j < (int)gp.N; ++j) {
			posPts[j].values[0] = x[j];// xFull[neigh_t + j*U];
			posPts[j].values[1] = y[j];// yFull[neigh_t + j*U];
			posPts[j].idx = j;
		}
		getNearestNeighbours(posPts, (int)gp.N, (int)GK, gp.L, fullGTENNidx);

		DataSet<GK + 2> dataForGTE((int)gp.N, (int)sp.KSG_neighbours);

		std::vector<double> neighbourAngles(GK);
		for(size_t i = 0; i < gp.N; i++) {
			dataForGTE.setData(h[i], (int)i, 0);
			for(size_t y_i = 0; y_i < (size_t)GK; ++y_i) {
				neighbourAngles[y_i] = h[fullGTENNidx[i*GK + y_i]];
			}
			//std::sort(neighbourAngles.begin(), neighbourAngles.end());
			for(size_t y_i = 0; y_i < (size_t)GK; ++y_i) {
				dataForGTE.setData(neighbourAngles[y_i], (int)i, (int)y_i + 1);
			}
			dataForGTE.setData(hnew[i], (int)i, GK + 1);
		}

		return calculateGTE(dataForGTE);
	}
	else {
		PEEXIT("Pairwise KSG not supported for ksg_gte_dims={}", sp.ksg_gte_dims);
	}
}

void fillAndRunPairwiseKSG(GeneratorParams const& gp,
						SimParams const& sp,
						std::vector<int> const& interactions,
						size_t const u,
						FlockState & fs,
						double const*const hAgg,
						int * fullGTENNidx,
						std::vector<VPoint<2> > & posPts,
						ITStats& itStats)
{
	auto const h = fs.h.reader(), x = fs.x.reader(), y = fs.y.reader();
	auto const hnew = fs.h.writer();
	if(sp.sims[MIBIN])
	{
		fmt::print("Calculating MI\n");
		double val = fillAndRunPairwiseMI(interactions, hnew, sp.KSG_neighbours);
		itStats.mi.accumulate(val);
	}

	if(sp.sims[TEBIN] /*&& u > 0*/ && u < sp.U - 1)
	{
		fmt::print("Calculating TE: ");
		double val = fillAndRunPairwiseTE(interactions, h, hnew, sp.KSG_neighbours);
		itStats.te.accumulate(val);
	}

	if(sp.sims[GTEBIN] /*&& u > 0*/ && u < sp.U - 1)
	{
		fmt::print("Calculating GTE:\n");
		// Get nearest neighbours for GTE to use
		double val = fillAndRunPairwiseGTE(gp, sp, posPts, x, y, h, hnew, hAgg, fullGTENNidx);
		itStats.gte.accumulate(val);
	}
}




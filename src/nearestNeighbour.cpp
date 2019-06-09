#include "../include/nearestNeighbour.hpp"
#include "../include/defs.hpp"
#include <fmt/format.h>
#include <queue>

void getNearestNeighbours(std::vector<VPoint<2> > &points, int const N, int const K, double const L, int * neighbours)
{
	int const K1 = K + 1;

	VpTreeEuclidean2D tree(L);
	tree.create(points);

	std::vector<VPoint<2> > results;
	std::vector<int> indices;
	std::vector<double> distances;

	for(int i = 0; i < N; ++i) {
		tree.knn_search(points[i], K1, &results, &indices, &distances);
		for(int k = 0; k < K; ++k) {
			neighbours[i*K + k] = results[k + 1].idx;
		}
	}
}

//Used for sanity checks
std::vector<int> brutedistkPoints(std::vector<VPoint<2> > const & dataPts, int const point, int const numPts, int const k, int const dim, double linearSize) {
	std::priority_queue<std::pair<double, int> > knn;
	for(int l = 0; l < k; ++l) {
		knn.push(std::pair<double, int>(linearSize, -1));
	}

	for(int j = 0; j < numPts; ++j) {
		double euclid = 0.0;
		//double maxnorm = 0.0;
		for(int d = 0; d < dim; ++d) {
			//maxnorm = std::max(maxnorm, dwrap(std::fabs(dataPts[point].values[d] - dataPts[j].values[d]), linearSize));
			double dist = dwrap(dataPts[point].values[d] - dataPts[j].values[d], linearSize);
			euclid += dist*dist;
		}
		euclid = sqrt(euclid);
		if(euclid < knn.top().first)
		//if(maxnorm < knn.top().first)
		{
			knn.push(std::pair<double, int>(euclid, j));
			//knn.push(std::pair<double, int>(maxnorm, j));
			if((signed int)knn.size() > k)
				knn.pop();
		}
	}
	std::vector<int> neighbours(k);
	int idx = k - 1;
	while(knn.empty() == false)
	{
		fmt::print("{}: {}\n", knn.top().second, knn.top().first);

		//neighbours.push_back(knn.top().second);
		neighbours[idx] = knn.top().second;
		idx--;
		knn.pop();
	}
	return neighbours;
}

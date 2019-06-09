#ifndef DATASET_BRUTEFORCE_HPP
#define DATASET_BRUTEFORCE_HPP

#include "vpUtil.hpp"
#include "psi.hpp"
#include "defs.hpp"
#include "vpTree.hpp"
#include "utils.hpp"
#include "dataSet.hpp"

#include <omp.h>
#include <random>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <fmt/format.h>

typedef double(*distFunc)(double, double);
double angleDist(double i, double j);

double regularDist(double i, double j);

/**
 * Holds data and functions for calculating dimensional quantites
 * needed by nearest-neighbour estimators (e.g., hyper-ball radii
 * and neighbour counting).
 * Templated on the number of dimensions stored in the dataset
 * Uses naive brute force searches for validating DataSet<D>
 */
template <unsigned int D>
class DataSetBruteForce
{
public:
	// Sets up data set holding N elements, using k nearest neighbours
	// for determining hyper-ball sizes
	DataSetBruteForce(int const N, int const k) : N(N), k(k)
	{
		ballsizes.resize(N);
		points.resize(N);
	}

	/**
	 * Copy settings and data from another dataset
	 */
	DataSetBruteForce(DataSet<D> const& other) :
		N(other.getN()), k(other.getK())
	{
		ballsizes.resize(N);
		points.resize(N);
		auto other_points = other.getPoints();
		std::copy(other_points, other_points + N, points.data());
	}

	// Get the elements being stored
	VPoint<D> const* getPoints() const {
		return points.data();
	}
	// Get the calculated hyper-ball sizes
	double const* getBallsizes() const {
		return ballsizes.data();
	}

	// Fill dimension d with data stream
	void setData(double const* const data, unsigned int const d) {
		if(d >= 0 && d < dim)
		{
			for(int i = 0; i < N; ++i) {
				points[i].values[d] = data[i];
				points[i].idx = i;
			}
		}
	}
	// Set element i of dimension d with data
	inline void setData(double const data, int i, unsigned int d) {
		points[i].values[d] = data;
		points[i].idx = i;
	}

	inline int getN() const { return N; }
	inline int getK() const { return k; }

	// Calculate the average ballsize of dataset
	double averageBallsize() {
		double total = std::accumulate(ballsizes.begin(), ballsizes.end(), 0.0,
			[](double const& acc, double const& x) {return acc + log(x * 2); });
		return (total/N);
	}

	// Calculate hyper-ball radii of dataset, using max-norm
	// distances
	void calcBallsizes() {
		if(k == 1) {
			for(int i = 0; i < N; ++i) {
				//int currentidx = 0;
				double currMaxNorm = DBL_MAX;
				for(int j = 0; j < N; ++j) {
					if(i == j) continue;
					double maxNormForJ = -1;
					for(unsigned int d = 0; d < dim; ++d) {
						double dist = dwrap(points[i].values[d] - points[j].values[d], TWOPI);
						maxNormForJ = std::max(maxNormForJ, dist);
					}
					if(maxNormForJ < currMaxNorm)
						currMaxNorm = maxNormForJ;
				}
				ballsizes[i] = currMaxNorm;
			}
		}
		else {
			for(int i = 0; i < N; ++i) {
				ballsizes[i] = brutedistk(points, i, N, k + 1, dim, angleDist);
			}
		}
	}

	// Calculate <psi(N_d(i)> for averaged every point along dimension D1
	// using the precalculated ballsizes. Optionally takes a pointer
	// to the ballsizes calculated via the non-naive method in DataSet
	double countNeighboursOneDim(const int d, double const*const ballsize_override = nullptr) {
		auto bs = ballsize_override ? ballsize_override : ballsizes.data();
		double total = 0.0;
		for(int i = 0; i < N; ++i) {
			int n = 0;
			double threshold = bs[i];
			for(int j = 0; j < N; ++j) {
				if(dwrap(points[i].values[d] - points[j].values[d], TWOPI) < threshold) {
					n++;
				}
			}
			if(n == 0)
				fmt::print("Counting N[1-bf] (i: {}, brute ballsizes: {}, ballsize: {}): n==0\n", i, ballsize_override ? "false" : "true", threshold);
			total += digamma(n);
		}
		return (total / N);
	}

	// Calculate <psi(N_d(i)> for every point along dimensions D1, D2
	// using the precalculated ballsizes. Optionally takes a pointer
	// to the ballsizes calculated via the non-naive method in DataSet
	double countNeighboursTwoDim(const int d1, const int d2, double const* const ballsize_override = nullptr) {
		auto bs = ballsize_override ? ballsize_override : ballsizes.data();
		double total = 0.0;
		for(int i = 0; i < N; ++i) {
			int n = 0;
			for(int j = 0; j < N; ++j) {
				double threshold = bs[i];
				if(dwrap(points[i].values[d1] - points[j].values[d1], TWOPI) < threshold &&
					dwrap(points[i].values[d2] - points[j].values[d2], TWOPI) < threshold) {
					n++;
				}
			}
			if(n == 0)
				fmt::print("Counting N[2-bf] (i={}, brute ballsizes : {}): n==0\n", i, ballsize_override ? "false" : "true");
			total += digamma(n);
		}
		return (total / N);
	}

	// Calculate <psi(N_d(i)> for every point along all dimensions except
	// the last using the precalculated ballsizes.
	// This is used for calculating N_xy, where y is multidimensional,
	// assuming dataset stores [X, Y', W]. Optionally takes a pointer
	// to the ballsizes calculated via the non-naive method in DataSet
	double countNeighboursNDim(const int d, double const* const ballsize_override = nullptr) {
		auto bs = ballsize_override ? ballsize_override : ballsizes.data();
		double total = 0.0;
		for(int i = 0; i < N; ++i) {
			int n = 0;
			for(int j = 0; j < N; ++j) {
				bool outside = false;
				for(int d1 = 0; d1 < d; ++d1) {
					double threshold = bs[i];
					if(dwrap(points[i].values[d1] - points[j].values[d1], TWOPI) >= threshold)
						outside = true;
				}
				if(outside == false)
					n++;
			}
			if(n == 0)
				fmt::print("Counting N[n-bf] (i={}, brute ballsizes : {}): n==0\n", i, ballsize_override ? "false" : "true");
			total += digamma(n);
		}
		return (total / N);
	}

	static unsigned int const dim = D;
private:
	std::vector<VPoint<D> > points;
	std::vector<double> ballsizes;
	int const N;
	int const k;
	
	// Finds the distance to the kth nearest neighbour using naive methods
	double brutedistk(std::vector<VPoint<D> >& dataPts, int const point, int const numPts, int const brute_k, int const brute_dim, distFunc dist) {
		std::priority_queue<std::pair<double, int> > knn;
		knn.push(std::pair<double, int>(DBL_MAX, -1));
		
		for(int j = 0; j < numPts; ++j) {
			double maxnorm = 0;
			for(int d = 0; d < brute_dim; ++d) {
				maxnorm = std::max(maxnorm, dist(dataPts[point].values[d], dataPts[j].values[d]));
			}
			if(maxnorm < knn.top().first) {
				knn.push(std::pair<double, int>(maxnorm, j));
				if((signed int)knn.size() > brute_k)
					knn.pop();
			}
		}
		return knn.top().first;
	}
};

// Compares the ballsizes calculated by both methods, and reports how
// many values have a difference greater than threshold.
// Will print out all differences if listDifferences flag is set.
template <unsigned int D>
bool compareBallsizeMethods(DataSet<D> const& dataSet, DataSetBruteForce<D> const& validationSet, double threshold, bool listDifferences = false) {
	bool same = true;
	auto const ballsizes = dataSet.getBallsizes();
	auto const validation_ballsizes = validationSet.getBallsizes();
	auto N = dataSet.getN();
	auto N1 = validationSet.getN();
	if (N != N1) {
		fmt::print("Cannot compare data sets. Differing number of elements: {} vs {}\n", N, N1);
		return false;
	}
	auto k = dataSet.getK(), k1 = validationSet.getK();
	if (k != k1) {
		fmt::print("Data set comparison meaningless. Differing k: {} vs {}\n", k, k1);
		return false;
	}
	for (int i = 0; i < N; ++i) {
		if (fabs(ballsizes[i] - validation_ballsizes[i]) > threshold) {
			if (listDifferences) {
				if (same) fmt::print("List of different ballsizes between optimised and brute force:\n"
					"id\topti\tbrute\n");
				fmt::print("{}\t{}\t{}\n", i, ballsizes[i], validation_ballsizes[i]);
			}
			same = false;
		}
	}
	return same;
}
#endif //DATASET_BRUTEFORCE_HPP

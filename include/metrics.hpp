#ifndef METRICS_HPP
#define METRICS_HPP

#include "psi.hpp"
#include "dataSet.hpp"

#include <fmt/format.h>

/**
 * Calculates MI, assuming X and Y are in dimensions 0 and 1 of
 * dataSet, respectively
 */
template <class Data>
double calculateMI(Data& dataSet) {
	static_assert(Data::dim >= 2, "Not enough dimensions to calculate MI");

	dataSet.calcBallsizes();
	double psiNX = dataSet.template countNeighboursOneDim<0>();
	double psiNY = dataSet.template countNeighboursOneDim<1>();

	return (psi(dataSet.getK()) + psi(dataSet.getN()) - psiNX - psiNY);
}

/**
 * Calculates TE, assuming X,Y,W are in dimensions 0,1,2 of
 * dataSet, respectively, where W = X+1
 */
template <class Data>
double calculateTE(Data& dataSet) {
	static_assert(Data::dim >= 3, "Not enough dimensions to calculate TE");

	dataSet.calcBallsizes();
	double psiNXW = dataSet.template countNeighboursTwoDim<0, 2>();
	double psiNXY = dataSet.template countNeighboursTwoDim<0, 1>();
	double psiNX = dataSet.template countNeighboursOneDim<0>();

	return (psi(dataSet.getK()) - psiNXW - psiNXY + psiNX);
}

/**
 * Calculate reduced GTE, using X and W, assuming X is in dim 0 and W is in dim 1
 * Where Tgl gets simplified to H(\Theta'_I|\Theta_I)-H(\Omega)
 */
template <class Data>
double calculateReducedGTE(Data& dataSet, double eta) {
	static_assert(Data::dim >= 2, "Not enough dimensions to calculate GTE");
	
	dataSet.calcBallsizes();
	double psiNX = dataSet.template countNeighboursOneDim<0>();
	double ballsizeAvg = dataSet.averageBallsize();

	return (-psi(dataSet.getK()) + psiNX + ballsizeAvg - log(eta));
}

/**
 * Calculates GTE for D dimensional dataSet. Assumes X is in
 * dimension 0 and W is in dimensions D-1, with Y consuming
 * dimensions 1 through to D-2
 */
template <class Data>
double calculateGTE(Data& dataSet) {
	static_assert(Data::dim >= 4, "Not enough dimensions to calculate GTE");

	dataSet.calcBallsizes();
	double psiNXW = dataSet.template countNeighboursTwoDim<0, Data::dim - 1>();
	double psiNXY = dataSet.countNeighboursNDim();
	double psiNX = dataSet.template countNeighboursOneDim<0>();

	return (psi(dataSet.getK()) - psiNXW - psiNXY + psiNX);
}


#endif //METRICS_HPP

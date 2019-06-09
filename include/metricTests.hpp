#ifndef METRICTESTS_HPP
#define METRICTESTS_HPP

#include "dataSet.hpp"
#include "dataSetBruteForce.hpp"
#include <fmt/format.h>

// Compares values calculated from each method and reports results
void compareValues(double opti, double brute_pure, double brute_opti,
				   char const*const metric, char const*const dim);

/**
 * Functions to verify accuracy of DataSet in calculating ballsizes and metrics.
 * Performs two verifications per metric, one against a dataset using completely
 * naive calculations, and a second using the ballsizes calculated by the DataSet
 * method rather than the naive calculation.
 */
template <class Data>
bool verifyMI(Data& dataSet, double threshold = 0.0) {
	static_assert(Data::dim >= 2, "Not enough dimensions to calculate MI");

	DataSetBruteForce<Data::dim> validation_dataset(dataSet);

	fmt::print("MI Verify: Starting verification\n");
	dataSet.calcBallsizes();
	auto const ballsizes = dataSet.getBallsizes();
	fmt::print("MI Verify: Ballsizes calculated\n");
	validation_dataset.calcBallsizes();
	fmt::print("MI Verify: Ballsizes calculated (validation)\n");

	bool ballsizesEqual = compareBallsizeMethods(dataSet, validation_dataset, threshold, true);
	fmt::print("MI Verify: Ballsizes {} equal\n", ballsizesEqual ? "are" : "are not");

	double psiNX = dataSet.template countNeighboursOneDim<0>();
	fmt::print("MI Verify: psi(n_x) calculated\n");
	double psiNXb = validation_dataset.countNeighboursOneDim(0);
	fmt::print("MI Verify: psi(n_x) calculated (validation)\n");
	double psiNXbb = validation_dataset.countNeighboursOneDim(0, ballsizes);
	fmt::print("MI Verify: psi(n_x) calculated (validation: using original ballsizes)\n");

	compareValues(psiNX, psiNXb, psiNXbb, "MI Verify", "n_x");

	double psiNY = dataSet.template countNeighboursOneDim<1>();
	fmt::print("MI Verify: psi(n_y) calculated\n");
	double psiNYb = validation_dataset.countNeighboursOneDim(1);
	fmt::print("MI Verify: psi(n_y) calculated (validation)\n");
	double psiNYbb = validation_dataset.countNeighboursOneDim(1, ballsizes);
	fmt::print("MI Verify: psi(n_y) calculated (validation: using original ballsizes)\n");

	compareValues(psiNY, psiNYb, psiNYbb, "MI Verify", "n_y");

	double opti = (psi(dataSet.getK()) + psi(dataSet.getN()) - psiNX - psiNY);
	double brute = (psi(dataSet.getK()) + psi(dataSet.getN()) - psiNXb - psiNYb);

	fmt::print("MI Verify:\nI : {}\nIb: {}\n", opti, brute);

	return (fabs(opti - brute) < threshold);
}

template <class Data>
bool verifyTE(Data& dataSet, double threshold = 0.0) {
	static_assert(Data::dim >= 3, "Not enough dimensions to calculate TE");

	DataSetBruteForce<Data::dim> validation_dataset(dataSet);

	dataSet.calcBallsizes();
	validation_dataset.calcBallsizes();

	bool ballsizesEqual = compareBallsizeMethods(dataSet, validation_dataset, threshold, true);
	if(!ballsizesEqual) {
		fmt::print("TE Verify: Ballsizes not equal\n");
	}
	auto const ballsizes = dataSet.getBallsizes();

	double psiNXW = dataSet.template countNeighboursTwoDim<0, 2>();
	double psiNXWb = validation_dataset.countNeighboursTwoDim(0, 2);
	double psiNXWbb = validation_dataset.countNeighboursTwoDim(0, 2, ballsizes);

	compareValues(psiNXW, psiNXWb, psiNXWbb, "TE Verify", "n_xw");

	double psiNXY = dataSet.template countNeighboursTwoDim<0, 1>();
	double psiNXYb = validation_dataset.countNeighboursTwoDim(0, 1);
	double psiNXYbb = validation_dataset.countNeighboursTwoDim(0, 1, ballsizes);

	compareValues(psiNXY, psiNXYb, psiNXYbb, "TE Verify", "n_xy");

	double psiNX = dataSet.template countNeighboursOneDim<0>();
	double psiNXb = validation_dataset.countNeighboursOneDim(0);
	double psiNXbb = validation_dataset.countNeighboursOneDim(0, ballsizes);

	compareValues(psiNX, psiNXb, psiNXbb, "TE Verify", "n_x");

	double opti = (psi(dataSet.getK()) - psiNXW - psiNXY + psiNX);
	double brute = (psi(dataSet.getK()) - psiNXWb - psiNXYb + psiNXb);

	fmt::print("TE Verify:\nI : {}\nIb: {}\n", opti, brute);
	return (fabs(opti - brute) < threshold);
}

template <class Data>
bool verifyGTE(Data& dataSet, double threshold = 0.0) {
	static_assert(Data::dim >= 4, "Not enough dimensions to calculate GTE");

	DataSetBruteForce<Data::dim> validation_dataset(dataSet);
	dataSet.calcBallsizes();
	validation_dataset.calcBallsizes();

	bool ballsizesEqual = compareBallsizeMethods(dataSet, validation_dataset, threshold, true);
	if(!ballsizesEqual) {
		fmt::print("GTE Verify: Ballsizes not equal\n");
	}
	auto const ballsizes = dataSet.getBallsizes();

	double psiNXW = dataSet.template countNeighboursTwoDim<0, Data::dim - 1>();
	double psiNXWb = validation_dataset.countNeighboursTwoDimBruteForce(0, Data::dim - 1);
	double psiNXWbb = validation_dataset.countNeighboursTwoDimBruteForce(0, Data::dim - 1, ballsizes);

	compareValues(psiNXW, psiNXWb, psiNXWbb, "GTE Verify", "n_xw");

	double psiNXY = dataSet.countNeighboursNDim();
	double psiNXYb = validation_dataset.countNeighboursNDimBruteForce(Data::dim - 1);
	double psiNXYbb = validation_dataset.countNeighboursNDimBruteForce(Data::dim - 1, ballsizes);

	compareValues(psiNXY, psiNXYb, psiNXYbb, "GTE Verify", "n_xy");

	//double psiNX = dataSet.countNeighboursOneDim(0);
	double psiNX = dataSet.template countNeighboursOneDim<0>();
	double psiNXb = validation_dataset.countNeighboursOneDimBruteForce(0);
	double psiNXbb = validation_dataset.countNeighboursOneDimBruteForce(0, ballsizes);

	compareValues(psiNX, psiNXb, psiNXbb, "GTE Verify", "n_x");

	double opti = (psi(dataSet.getK()) - psiNXW - psiNXY + psiNX);
	double brute = (psi(dataSet.getK()) - psiNXWb - psiNXYb + psiNXb);

	fmt::print("GTE Verify:\nI : {}\nIb: {}\n", opti, brute);
	return (fabs(opti - brute) < threshold);
}

#endif //METRICTESTS_HPP
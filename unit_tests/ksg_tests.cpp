#include <catch2/catch.hpp>

// Already have test functions in dataSet.hpp and metricTests.hpp
// May change later to use unit tests instead

//Write tests that take advantage of the tests already embedded in metric tests...derp

#include "../include/dataSet.hpp"
#include "../include/metricTests.hpp"
#include "../include/metrics.hpp"
#include <random>

TEST_CASE("Dataset calculates MI correctly", "[KSG]") {
	std::mt19937_64 rng(1);
	std::uniform_real_distribution<double> dis;

	SECTION("N=1000, K_K=3") {
		const int N = 1000;
		const int K = 3;
		DataSet<2> data(N, K);

		for(int i = 0; i < N; ++i) {
			data.setData(PI*(2.0*dis(rng) - 1.0), i, 0);
			data.setData(PI*(2.0*dis(rng) - 1.0), i, 1);
		}

		REQUIRE(verifyMI(data, 0.000001));
	}

	SECTION("Shifting by PI/2") {
		const int N = 1000;
		const int K = 3;
		DataSet<2> data(N, K);
		DataSet<2> dataShifted(N, K);

		for(int i = 0; i < N; ++i) {
			double x = PI*(2.0*dis(rng) - 1.0);
			double y = PI*(2.0*dis(rng) - 1.0);
			data.setData(x, i, 0);
			data.setData(y, i, 1);
			dataShifted.setData(awrap(x + PI / 2), i, 0);
			dataShifted.setData(awrap(y + PI / 2), i, 1);
		}

		REQUIRE(calculateMI(data) == Approx(calculateMI(dataShifted)));
	}

}

TEST_CASE("Dataset calculates TE/Avg GTE correctly", "[KSG]") {
	std::mt19937_64 rng(1);
	std::uniform_real_distribution<double> dis;

	SECTION("N=1000, K_K=3") {
		const int N = 1000;
		const int K = 3;
		DataSet<3> data(N, K);

		double prev = PI*(2.0*dis(rng) - 1.0);
		for(int i = 0; i < N; ++i) {
			double x = PI*(2.0*dis(rng) - 1.0);
			data.setData(x, i, 0);
			data.setData(PI*(2.0*dis(rng) - 1.0), i, 1);
			data.setData(prev, i, 2);
			prev = x;
		}

		REQUIRE(verifyTE(data, 0.000001));
	}

}
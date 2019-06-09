#include <catch2/catch.hpp>

#include "../include/model.hpp"
#include "../include/vpUtil.hpp"
#include "../include/vpTree.hpp"
#include "../include/utils.hpp"
//#include "../src/model.cpp"

TEST_CASE("Vicsek metric simulates boids", "[vicsek]") {
	std::mt19937_64 rng(1);
	std::uniform_real_distribution<double> dis;

	SECTION("Single particle case") {
		size_t N = 1;
		FlockState fs(N);
		double L = 10;
		double v = 0.1;
		fs.h.writer()[0] = 0;
		fs.x.writer()[0] = 5;
		fs.y.writer()[0] = 5;
		fs.swapBuffers();
		double hAgg;
		VpTreeEuclidean2D vicsekTree(L);

		SECTION("No noise = no angle change and particle moves correct distance") {
			vicsekMetric(N, L, v, 0.0, fs, &hAgg, rng, 0, vicsekTree);
			REQUIRE(fs.h.reader()[0] == Approx(fs.h.writer()[0]));
			REQUIRE(dsqwrap(fs.x.reader()[0] - fs.x.writer()[0], fs.y.reader()[0] - fs.y.writer()[0], L) == Approx(v*v));
			REQUIRE(hAgg == Approx(0));
		}
		SECTION("Rng gets reset in subsections") {
			INFO("First number in random number stream with seed = 1 is 0.1338766");
			REQUIRE(dis(rng) == Approx(0.1338766));
		}
		SECTION("2pi noise = completely new angle and particle moves correct distance") {
			vicsekMetric(N, L, v, TWOPI, fs, &hAgg, rng, 0, vicsekTree);
			REQUIRE(fs.h.writer()[0] >= -PI);
			REQUIRE(fs.h.writer()[0] < PI);
			REQUIRE(dsqwrap(fs.x.reader()[0] - fs.x.writer()[0], fs.y.reader()[0] - fs.y.writer()[0], L) == Approx(v*v));
			REQUIRE(hAgg == Approx(0));
		}
	}

	SECTION("Two particle case - Nearby") {
		const size_t N = 2;
		double L = 10;
		double v = 0.1;
		double e = 0.1;
		FlockState fs(N);
		fs.h.writer()[0] = -e;
		fs.h.writer()[1] = e;
		fs.x.writer()[0] = 5;
		fs.x.writer()[1] = 5.1;
		fs.y.writer()[0] = 5;
		fs.y.writer()[1] = 5.1;
		fs.swapBuffers();
		double hAgg[2];
		VpTreeEuclidean2D vicsekTree(L);

		SECTION("No noise = new angle for both particles should cancel to 0") {
			INFO("Particle 0 has angle " << -e << " while particle 1 has angle " << e << " which when combined should be 0");
			vicsekMetric(N, L, v, 0.0, fs, hAgg, rng, 0, vicsekTree);
			REQUIRE(fs.h.writer()[0] == Approx(fs.h.writer()[1]));
			REQUIRE(fs.h.writer()[0] == Approx(0));
			REQUIRE(dsqwrap(fs.x.reader()[0] - fs.x.writer()[0], fs.y.reader()[0] - fs.y.writer()[0], L) == Approx(v*v));
			REQUIRE(dsqwrap(fs.x.reader()[1] - fs.x.writer()[1], fs.y.reader()[1] - fs.y.writer()[1], L) == Approx(v*v));
			REQUIRE(hAgg[0] == Approx(0));
			REQUIRE(hAgg[1] == Approx(0));
		}

		SECTION("Bordering boundary") {
			FlockState fs2(N);
			fs2.h.writer()[0] = -e;
			fs2.h.writer()[1] = e;
			fs2.x.writer()[0] = 9.95;
			fs2.x.writer()[1] = 0.05;
			fs2.y.writer()[0] = 9.95;
			fs2.y.writer()[1] = 0.05;
			fs2.swapBuffers();
			double hAgg2[2];
			vicsekMetric(N, L, v, 0.0, fs, hAgg, rng, 0, vicsekTree);
			vicsekMetric(N, L, v, 0.0, fs2, hAgg2, rng, 0, vicsekTree);

			REQUIRE(fs.h.writer()[0] == Approx(fs2.h.writer()[0]));
			REQUIRE(fs.h.writer()[1] == Approx(fs2.h.writer()[1]));
			REQUIRE(hAgg[0] == Approx(hAgg2[0]));
			REQUIRE(hAgg[1] == Approx(hAgg2[1]));
			REQUIRE(dsqwrap(fs2.x.reader()[0] - fs2.x.writer()[0], fs2.y.reader()[0] - fs2.y.writer()[0], L) == Approx(v*v));
			REQUIRE(dsqwrap(fs2.x.reader()[1] - fs2.x.writer()[1], fs2.y.reader()[1] - fs2.y.writer()[1], L) == Approx(v*v));

		}
	}

	SECTION("Two particle case - Far") {
		const size_t N = 2;
		FlockState fs(N);
		double L = 10;
		double v = 0.1;
		double e = 0.1;
		fs.h.writer()[0] = -e;
		fs.h.writer()[1] = e;
		fs.x.writer()[0] = 3;
		fs.x.writer()[1] = 6;
		fs.y.writer()[0] = 3;
		fs.y.writer()[1] = 6;
		fs.swapBuffers();
		double hAgg[2];
		VpTreeEuclidean2D vicsekTree(L);

		SECTION("No noise = new angle for both particles should be their old angles") {
			INFO("Particle 0 has angle " << -e << " while particle 1 has angle " << e << " which won't combine due to interaction radius");
			vicsekMetric(N, L, v, 0.0, fs, hAgg, rng, 0, vicsekTree);
			REQUIRE(fs.h.writer()[0] == Approx(fs.h.reader()[0]));
			REQUIRE(fs.h.writer()[1] == Approx(fs.h.reader()[1]));
			REQUIRE(dsqwrap(fs.x.reader()[0] - fs.x.writer()[0], fs.y.reader()[0] - fs.y.writer()[0], L) == Approx(v*v));
			REQUIRE(dsqwrap(fs.x.reader()[1] - fs.x.writer()[1], fs.y.reader()[1] - fs.y.writer()[1], L) == Approx(v*v));
			REQUIRE(hAgg[0] == Approx(fs.h.reader()[0]));
			REQUIRE(hAgg[1] == Approx(fs.h.reader()[1]));
		}
	}

	SECTION("Many particle case is true to brute force version of model") {
		const size_t N = 100;
		double L = 10;
		double v = 0.1;
		FlockState fs(N);
		double* s = new double[N];
		double* h = new double[N];
		double* hAggOpti = new double[N];
		double* x = new double[N];
		double* y = new double[N];
		double* hnewBrute = new double[N];
		double* xnewBrute = new double[N];
		double* ynewBrute = new double[N];
		VpTreeEuclidean2D vicsekTree(L);

		initialise(0, N, L, h, x, y, s, rng);
		auto h1 = fs.h.writer(), x1 = fs.x.writer(), y1 = fs.y.writer();
		for (size_t i = 0; i < N; ++i) {
			h1[i] = h[i];
			x1[i] = x[i];
			y1[i] = y[i];
		}
		fs.swapBuffers();
		rng.seed(1);
		std::mt19937_64 rngBrute(1);

		REQUIRE(dis(rng) == Approx(dis(rngBrute)));

		SECTION("No noise simulation is accurate to brute force") {
			const double eta = 0.0;
			vicsek(N, L, v, eta, h, x, y, hnewBrute, xnewBrute, ynewBrute,
				   rngBrute);
			vicsekMetric(N, L, v, eta, fs, hAggOpti, rng, 0, vicsekTree);

			for(size_t i = 0; i < N; ++i) {
				CAPTURE(i);
				REQUIRE(fs.h.writer()[i] == Approx(hnewBrute[i]));
				REQUIRE(fs.x.writer()[i] == Approx(xnewBrute[i]));
				REQUIRE(fs.y.writer()[i] == Approx(ynewBrute[i]));
			}
		}

		INFO("Can't test noisy variants as the optimisation causes the particles " \
			 "to be processed in a different order, and thus get a different random number stream");

		delete[] h;
		delete[] x;
		delete[] y;
		delete[] hAggOpti;
		delete[] xnewBrute;
		delete[] ynewBrute;
		delete[] hnewBrute;
	}
}

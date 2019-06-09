#include <catch2/catch.hpp>

#include "../include/defs.hpp"

TEST_CASE("Wrapping functions wrap properly", "[utils]") {
	SECTION("Coordinate wrapping between [0,L) works properly") {
		REQUIRE(cwrap(25, 50) == Approx(25));
		REQUIRE(cwrap(75, 50) == Approx(25));
		REQUIRE(cwrap(-25, 50) == Approx(25));
	}	

	SECTION("Coordinate wrapping between [-L/2, L/2) works properly") {
		REQUIRE(rwrap(10, 50) == Approx(10));
		REQUIRE(rwrap(60, 50) == Approx(10));
		REQUIRE(rwrap(25, 50) == Approx(-25));
		REQUIRE(rwrap(30, 50) == Approx(-20));
		REQUIRE(rwrap(50, 50) == Approx(0));
	}

	SECTION("Getting the distance in linear space gets shortest distance") {
		REQUIRE(dwrap(2, 10) == Approx(2));
		REQUIRE(dwrap(8, 10) == Approx(2));
		REQUIRE(dwrap(-2, 10) == Approx(2));
		REQUIRE(dwrap(-8, 10) == Approx(2));
	}

	SECTION("Angle gets wrapped between [-pi, pi) properly") {
		REQUIRE(awrap(0) == Approx(0));
		REQUIRE(awrap(PI) == Approx(-PI));
		REQUIRE(awrap(-PI) == Approx(-PI));
		double e = 0.2;
		REQUIRE(awrap(PI + e) == Approx(-PI + e));
		REQUIRE(awrap(-PI - e) == Approx(PI - e));
		REQUIRE(awrap(e) == Approx(e));
		REQUIRE(awrap(-e) == Approx(-e));
	}

	SECTION("Angle gets binned properly") {
		REQUIRE(binit(-PI, 16) == 0);
		REQUIRE(binit(0, 16) == 8);
		REQUIRE(binit(0.5 * (TWOPI / 16.0), 16) == 8);
		REQUIRE(binit(PI, 16) == 15);
		REQUIRE(binit(-PI, 64) == 0);
		REQUIRE(binit(0, 64) == 32);
		REQUIRE(binit(0.5 * (TWOPI / 64.0), 64) == 32);
		REQUIRE(binit(-PI, 10) == 0);
		REQUIRE(binit(0.1, 10) == 5);
	}
	
	SECTION("Angle gets binned properly (binit2)") {
		REQUIRE(binit2(-TWOPI, 16) == 0);
		REQUIRE(binit2(0, 16) == 8);
		REQUIRE(binit2(0.5 * (2*TWOPI / 16.0), 16) == 8);
		REQUIRE(binit2(TWOPI, 16) == 15);
		REQUIRE(binit2(-TWOPI, 64) == 0);
		REQUIRE(binit2(0, 64) == 32);
		REQUIRE(binit2(0.5 * (2*TWOPI / 64.0), 64) == 32);
		REQUIRE(binit2(-TWOPI, 10) == 0);
		REQUIRE(binit2(0.1, 10) == 5);
	}
	
	SECTION("Angle gets binned properly (binit3") {
		REQUIRE(binit3(0, 16) == 0);
		REQUIRE(binit3(0.5 * PI, 16) == 8);
		REQUIRE(binit3(PI, 16) == 15);
		REQUIRE(binit3(0, 64) == 0);
		REQUIRE(binit3(0.5 * PI, 64) == 32);
		REQUIRE(binit3(PI, 10) == 9);
	}
}

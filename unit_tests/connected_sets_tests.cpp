#include <catch2/catch.hpp>

#include "../include/connectedFlocks.hpp"

TEST_CASE("Connected sets algorithm works properly", "[utils]") {
	SECTION("Disjoint sets") {
		ConnectedSets connectedSet(10);

		for(int i = 0; i < 10; ++i)
			connectedSet.add(i);

		for(int i = 0; i < 10; ++i) {
			REQUIRE(connectedSet.particleMap[i] == i);
			REQUIRE(connectedSet.sets[i].size() == 1);
		}
	}

	SECTION("Duplicates") {
		ConnectedSets connectedSet(10);

		connectedSet.connect(0, 1);
		connectedSet.connect(0, 1);
		connectedSet.connect(0, 1);
		connectedSet.connect(1, 0);
		connectedSet.connect(1, 0);
		connectedSet.connect(1, 0);
		connectedSet.connect(0, 1);
		connectedSet.connect(1, 0);

		REQUIRE(connectedSet.sets[0].size() == 2);
		REQUIRE(connectedSet.sets[1].size() == 0);
	}

	SECTION("Duplicates - Merging") {
		ConnectedSets connectedSet(10);

		connectedSet.connect(0, 1);
		connectedSet.connect(0, 1);
		connectedSet.connect(2, 3);
		connectedSet.connect(2, 3);
		connectedSet.connect(2, 3);
		connectedSet.connect(1, 0);
		connectedSet.connect(1, 0);
		connectedSet.connect(1, 0);
		connectedSet.connect(0, 1);
		connectedSet.connect(3, 2);
		connectedSet.connect(1, 0);
		connectedSet.connect(1, 2);
		connectedSet.connect(3, 0);

		REQUIRE(connectedSet.sets[0].size() == 4);
		REQUIRE(connectedSet.sets[1].size() == 0);
	}

	SECTION("Basic usage - no merging") {
		ConnectedSets connectedSet(10);

		connectedSet.connect(0, 1);
		connectedSet.connect(2, 3);
		connectedSet.connect(1, 4);
		connectedSet.connect(3, 5);
		connectedSet.connect(0, 6);
		connectedSet.connect(6, 7);

		REQUIRE(connectedSet.sets[0].size() == 5);
		REQUIRE(connectedSet.sets[1].size() == 3);

		REQUIRE(connectedSet.particleMap[0] == 0);
		REQUIRE(connectedSet.particleMap[1] == 0);
		REQUIRE(connectedSet.particleMap[2] == 1);
		REQUIRE(connectedSet.particleMap[3] == 1);
		REQUIRE(connectedSet.particleMap[4] == 0);
		REQUIRE(connectedSet.particleMap[5] == 1);
		REQUIRE(connectedSet.particleMap[6] == 0);
		REQUIRE(connectedSet.particleMap[7] == 0);
	}

	SECTION("Basic usage - merging") {
		ConnectedSets connectedSet(10);

		connectedSet.connect(0, 1);
		connectedSet.connect(2, 3);
		connectedSet.connect(4, 5);
		connectedSet.connect(1, 2);
		connectedSet.connect(2, 6);

		REQUIRE(connectedSet.sets[0].size() == 5);
		REQUIRE(connectedSet.sets[1].size() == 0);
		REQUIRE(connectedSet.sets[2].size() == 2);

		REQUIRE(connectedSet.particleMap[0] == 0);
		REQUIRE(connectedSet.particleMap[1] == 0);
		REQUIRE(connectedSet.particleMap[2] == 0);
		REQUIRE(connectedSet.particleMap[3] == 0);
		REQUIRE(connectedSet.particleMap[4] == 2);
		REQUIRE(connectedSet.particleMap[5] == 2);
		REQUIRE(connectedSet.particleMap[6] == 0);
		REQUIRE(connectedSet.holes.size() == 1);
		REQUIRE(connectedSet.holes.top() == 1);

		connectedSet.squash();
		REQUIRE(connectedSet.holes.empty());
		REQUIRE(connectedSet.particleMap[4] == 1);
		REQUIRE(connectedSet.particleMap[5] == 1);
		REQUIRE(connectedSet.sets[1].size() == 2);
		REQUIRE(connectedSet.sets[2].size() == 0);

		SECTION("Advanced usage - clear") {
			connectedSet.clear();

			REQUIRE(connectedSet.nextSet == 0);

			connectedSet.connect(0, 1);
			connectedSet.connect(2, 3);
			connectedSet.connect(1, 4);
			connectedSet.connect(3, 5);
			connectedSet.connect(0, 6);
			connectedSet.connect(6, 7);

			REQUIRE(connectedSet.sets[0].size() == 5);
			REQUIRE(connectedSet.sets[1].size() == 3);

			REQUIRE(connectedSet.particleMap[0] == 0);
			REQUIRE(connectedSet.particleMap[1] == 0);
			REQUIRE(connectedSet.particleMap[2] == 1);
			REQUIRE(connectedSet.particleMap[3] == 1);
			REQUIRE(connectedSet.particleMap[4] == 0);
			REQUIRE(connectedSet.particleMap[5] == 1);
			REQUIRE(connectedSet.particleMap[6] == 0);
			REQUIRE(connectedSet.particleMap[7] == 0);
		}
	}

	SECTION("Basic usage - large squash") {
		ConnectedSets connectedSet(14);
		connectedSet.connect(0, 1);
		connectedSet.connect(2, 3);
		connectedSet.connect(4, 5);
		connectedSet.connect(6, 7);
		connectedSet.connect(8, 9);
		connectedSet.connect(10, 11);
		connectedSet.connect(12, 13);

		connectedSet.connect(1, 2);
		connectedSet.connect(5, 6);
		connectedSet.connect(9, 10);

		REQUIRE(connectedSet.sets[1].size() == 0);
		REQUIRE(connectedSet.sets[3].size() == 0);
		REQUIRE(connectedSet.sets[5].size() == 0);

		connectedSet.squash();

		REQUIRE(connectedSet.holes.empty());
		REQUIRE(connectedSet.particleMap[0] == 0);
		REQUIRE(connectedSet.particleMap[1] == 0);
		REQUIRE(connectedSet.particleMap[2] == 0);
		REQUIRE(connectedSet.particleMap[3] == 0);
		REQUIRE(connectedSet.particleMap[4] == 2);
		REQUIRE(connectedSet.particleMap[5] == 2);
		REQUIRE(connectedSet.particleMap[6] == 2);
		REQUIRE(connectedSet.particleMap[7] == 2);
		REQUIRE(connectedSet.particleMap[8] == 3);
		REQUIRE(connectedSet.particleMap[9] == 3);
		REQUIRE(connectedSet.particleMap[10] == 3);
		REQUIRE(connectedSet.particleMap[11] == 3);
		REQUIRE(connectedSet.particleMap[12] == 1);
		REQUIRE(connectedSet.particleMap[13] == 1);
	}

	SECTION("Basic usage - large squash (last set is hole)") {
		ConnectedSets connectedSet(16);
		connectedSet.connect(0, 1);
		connectedSet.connect(2, 3);
		connectedSet.connect(4, 5);
		connectedSet.connect(6, 7);
		connectedSet.connect(8, 9);
		connectedSet.connect(10, 11);
		connectedSet.connect(12, 13);
		connectedSet.connect(14, 15);

		connectedSet.connect(1, 2);
		connectedSet.connect(5, 6);
		connectedSet.connect(9, 10);
		connectedSet.connect(13, 14);

		REQUIRE(connectedSet.sets[1].size() == 0);
		REQUIRE(connectedSet.sets[3].size() == 0);
		REQUIRE(connectedSet.sets[5].size() == 0);
		REQUIRE(connectedSet.sets[7].size() == 0);

		connectedSet.squash();

		REQUIRE(connectedSet.holes.empty());
		/*REQUIRE(connectedSet.particleMap[0] == 0);
		REQUIRE(connectedSet.particleMap[1] == 0);
		REQUIRE(connectedSet.particleMap[2] == 0);
		REQUIRE(connectedSet.particleMap[3] == 0);
		REQUIRE(connectedSet.particleMap[4] == 2);
		REQUIRE(connectedSet.particleMap[5] == 2);
		REQUIRE(connectedSet.particleMap[6] == 2);
		REQUIRE(connectedSet.particleMap[7] == 2);
		REQUIRE(connectedSet.particleMap[8] == 3);
		REQUIRE(connectedSet.particleMap[9] == 3);
		REQUIRE(connectedSet.particleMap[10] == 3);
		REQUIRE(connectedSet.particleMap[11] == 3);
		REQUIRE(connectedSet.particleMap[12] == 1);
		REQUIRE(connectedSet.particleMap[13] == 1);
		REQUIRE(connectedSet.particleMap[14] == 1);
		REQUIRE(connectedSet.particleMap[15] == 1);*/
	}
}
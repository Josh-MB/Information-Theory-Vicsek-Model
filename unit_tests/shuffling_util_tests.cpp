#include <catch2/catch.hpp>

#include "../include/defs.hpp"
#include <random>

TEST_CASE("Shuffling algorithm shuffles correctly", "[utils]") {
	std::mt19937_64 rng(1);

	SECTION("Only shuffles, doesn't delete/duplicate items") {
		std::vector<InteractionPair> indices;
		int check[100];

		for(int i = 0; i < 100; ++i) {
			indices.push_back(InteractionPair(i, i));
			check[i] = 0;
		}

		shuffleFirstMIndices(indices, 10, rng);

		for(int i = 0; i < 100; ++i) {
			check[indices[i].first]++;
		}

		for(int i = 0; i < 100; ++i) {
			REQUIRE(check[i] == 1);
		}

	}
}
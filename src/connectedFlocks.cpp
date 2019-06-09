#include "../include/connectedFlocks.hpp"
#include "../include/stats.hpp"
#include "../include/utils.hpp"

#include <fmt/format.h>
#include <cmath>

void calcAverageFlockDirections(ConnectedSets& connectedSets, double const* h, std::vector<double>& avg)
{
	if(h == 0)
		return;

	avg.clear();

	for(size_t i = 0; i < connectedSets.sets.size(); ++i) {
		double xavg = 0.0, yavg = 0.0;

		for(size_t j = 0; j < connectedSets.sets[i].size(); ++j) {
			double x, y;
			sincos(h[connectedSets.sets[i][j]], &y, &x);
			xavg += x;
			yavg += y;
		}
		avg.push_back(awrap(atan2(yavg, xavg)));
	}
}

void calcFlockCenters(ConnectedSets& connectedSets, double const*const x, double const*const y, std::vector<double>& cm_x, std::vector<double>& cm_y, const double L)
{
	if(x == 0 || y == 0)
		return;

	cm_x.clear();
	cm_y.clear();

	for(size_t set = 0; set < connectedSets.sets.size(); ++set) {
		std::vector<double> temp_x, temp_y;

		for(size_t i = 0; i < connectedSets.sets[set].size(); ++i) {
			temp_x.push_back(x[connectedSets.sets[set][i]]);
			temp_y.push_back(y[connectedSets.sets[set][i]]);
		}

		double cx, cy;
		calcCenterMass(connectedSets.sets[set].size(), &temp_x[0], &temp_y[0], L, cx, cy);
		cm_x.push_back(cx);
		cm_y.push_back(cy);
	}
}

void calcFlockRadii(ConnectedSets & connectedSets, double const * const x, double const * const y, std::vector<double> const & cm_x, std::vector<double> const & cm_y, const double L, std::vector<double>& radii)
{
	if(x == 0 || y == 0)
		return;

	radii.clear();

	for(size_t set = 0; set < connectedSets.sets.size(); ++set) {
		double best_mag = 0;

		for(size_t i = 0; i < connectedSets.sets[set].size(); ++i) {
			double temp_x = offsetwrap(cm_x[set], x[connectedSets.sets[set][i]], L);
			double temp_y = offsetwrap(cm_y[set], y[connectedSets.sets[set][i]], L);
			double mag = temp_x*temp_x + temp_y*temp_y;
			best_mag = std::max(mag, best_mag);
		}

		radii.push_back(std::sqrt(best_mag));
	}
}

void printFlock(ConnectedSets& connectedSet)
{
	fmt::print("Connected flocks\n");
	for(size_t i = 0; i < connectedSet.sets.size(); ++i) {
		//if(connectedSet.sets[i].size() != 0) {
			fmt::print("Flock {} (size: {}):\n", i, connectedSet.sets[i].size());
			for(size_t j = 0; j < connectedSet.sets[i].size(); ++j) {
				fmt::print("{}, ", connectedSet.sets[i][j]);
			}
			fmt::print("\n");
		//}
	}

	fmt::print("\nParticle map\n");
	for(size_t i = 0; i < connectedSet.particleMap.size(); ++i) {
		fmt::print("{} -> {}\n", i, connectedSet.particleMap[i]);
	}
	fmt::print("\n\n");
}

int getLargestFlock(ConnectedSets& connectedSets, int& flockSize)
{
	size_t largest = 0;
	size_t idx = 0;
	for(size_t i = 0; i < connectedSets.sets.size(); ++i) {
		if(connectedSets.sets[i].size() > largest) {
			largest = connectedSets.sets[i].size();
			idx = i;
		}
	}
	flockSize = (int)largest;
	return (int)idx;
}

void serialize(ConnectedSets& connectedSet, FILE* fp)
{
	size_t numSets = connectedSet.sets.size();
	size_t numParticles = connectedSet.particleMap.size();
	fwrite(&numSets, sizeof(size_t), 1, fp);
	fwrite(&numParticles, sizeof(size_t), 1, fp);
	fwrite(&(connectedSet.particleMap[0]), sizeof(int), numParticles, fp);
	for(size_t i = 0; i < numSets; ++i) {
		size_t setSize = connectedSet.sets[i].size();
		fwrite(&setSize, sizeof(size_t), 1, fp);
		fwrite(&(connectedSet.sets[i][0]), sizeof(int), setSize, fp);
	}
}

void load(ConnectedSets& connectedSet, FILE* fp)
{
	size_t numSets, numParticles;
	if(fread(&numSets, sizeof(size_t), 1, fp) != 1) PEEXIT("read failed");
	if(fread(&numParticles, sizeof(size_t), 1, fp) != 1) PEEXIT("read failed");
	connectedSet.particleMap.resize(numParticles);
	if(fread(&(connectedSet.particleMap[0]), sizeof(int), numParticles, fp) != numParticles) PEEXIT("read failed");
	connectedSet.sets.resize(numSets);
	for(size_t i = 0; i < numSets; ++i) {
		connectedSet.sets[i].reserve(numParticles);
		size_t setSize;
		if(fread(&setSize, sizeof(size_t), 1, fp) != 1) PEEXIT("read failed");
		connectedSet.sets[i].resize(setSize);
		if(fread(&(connectedSet.sets[i][0]), sizeof(int), setSize, fp) != setSize) PEEXIT("read failed");
	}
}

#include "defs.hpp"
#include "utils.hpp"

#include <vector>

/**
 * Manage histograms for calculating information theory
 * metrics from discretised angles. Calculation not done
 * here. Histogram written to file for external processing.
 */
class ITHistograms
{
private:
	std::vector<uint64_t> miHist;
	std::vector<uint64_t> marginalMiHist;
	std::vector<uint64_t> teHist;
	std::vector<uint64_t> gteHist;
	std::vector<size_t> binBuffer;
	size_t B;
	int gte_dims;
	bool sims[SIM_COUNT] = { false };
	bool shuffle;
public:
	ITHistograms(GeneratorParams const& gp, SimParams const& sp);

	void accumulateHists(GeneratorParams const& gp,
		FlockState & fs,
		double const rotation,
		double const cumulRotation,
		std::vector<int> const & interactionSet
	);

	void writeToFile(FILE *fp, SIM sim);
};
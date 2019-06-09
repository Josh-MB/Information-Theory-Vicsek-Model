#include "../include/itHistograms.hpp"
#include "../include/stats.hpp"

ITHistograms::ITHistograms(GeneratorParams const& gp, SimParams const& sp) :
	B(sp.B), gte_dims(sp.hist_gte_dims), shuffle(sp.shuffle)
{
	for (int i = 0; i < SIM_COUNT; ++i) {
		sims[i] = sp.sims[i];
	}
	
	if (sp.sims[MIBIN]) {
		miHist.resize(B);
		marginalMiHist.resize(B);
	}
	if (sp.sims[TEBIN])
		teHist.resize(B*B*B);
	if (sp.sims[GTEBIN]) {
		if (sp.hist_gte_dims == 1) {
			gteHist.resize(B);
		}
		else if (sp.hist_gte_dims == 2) {
			gteHist.resize(B * B);
		}
	}
	binBuffer.resize(gp.N);
}

void ITHistograms::accumulateHists(GeneratorParams const & gp, FlockState & fs, double const rotation, double const cumulRotation, std::vector<int> const & interactionSet)
{
	auto h = fs.h.reader(), x = fs.x.reader(), y = fs.y.reader();
	auto hnew = fs.h.writer();
	if (sims[MIBIN] && !shuffle) {
		//MI_histogram(p.N, h, cumulativeRotations[u], p.B, miHist, bins, allInteractions[u]);
		MI_1D_histogram(gp.N, gp.L, h, x, y, B, miHist.data(), marginalMiHist.data());
	}

	if (sims[TEBIN] && !shuffle)
		TE_histogram(gp.N, h, hnew, rotation, cumulRotation, B, teHist.data(), binBuffer.data(), interactionSet);
	if (sims[GTEBIN] && !shuffle)// && u > 1)
	{
		if (gte_dims == 1) {
			//GTE_avg_histogram(p.N, h, hnew, rotations[u], cumulativeRotations[u], hAgg, p.B, gteHist);
			GTE_1D_histogram(gp.N, h, hnew, B, gteHist.data(), rotation);
		}
		else if (gte_dims == 2) {
			GTE_2D_histogram(gp.N, h, hnew, B, gteHist.data(), rotation);
		}
	}
}

void ITHistograms::writeToFile(FILE * fp, SIM sim)
{
	switch (sim) {
	case MIBIN:
		//if(fwrite(miHist, sizeof(uint64_t), p.B*p.B, fp) != p.B*p.B) PEEXIT("writing MI histogram failed");
		if (fwrite(miHist.data(), sizeof(uint64_t), B, fp) != B) PEEXIT("writing MI histogram failed");
		if (fwrite(marginalMiHist.data(), sizeof(uint64_t), B, fp) != B) PEEXIT("writing MI marginal histogram failed");
		break;
	case TEBIN:
		if (fwrite(teHist.data(), sizeof(uint64_t), B*B*B, fp) != B*B*B) PEEXIT("writing TE histogram failed");
		break;
	case GTEBIN:
		if (gte_dims == 1) {
			//if(fwrite(gteHist, sizeof(uint64_t), B*B*B, fp) != B*B*B) PEEXIT("writing GTE histogram failed");
			if (fwrite(gteHist.data(), sizeof(uint64_t), B, fp) != B) PEEXIT("writing GTE histogram failed");
		}
		else if (gte_dims == 2) {
			if (fwrite(gteHist.data(), sizeof(uint64_t), B*B, fp) != B*B) PEEXIT("writing 2D GTE histogram failed");
		}
		break;
	default:
		break;
	}
}
//
// Created by max on 6/27/16.
//

#ifndef RNARK_AALBERTS_SCORER_HPP
#define RNARK_AALBERTS_SCORER_HPP

#include "nn_scorer.hpp"
#include "multi_array.hpp"
#include "aalberts_model.hpp"
namespace librnary {


class AalbertsScorer: public NNScorer<AalbertsModel> {
private:
	/**
	 * Returns a DP table containing the optimal configurations of stacking and a/b length segments.
	 * dp_table[as][bs][s] is the MFE configuration from subsurface s onward.
	 * The configuration contains as  a-length segments, and bs b-length segments, in the prefix preceding s.
	 * Optimal "stacking" is defined here http://rna.urmc.rochester.edu/NNDB/turner04/exterior.html
	 * Aalberts & Nandagopal (2010) define a and b length segments.
	 * Takes a list of features (branches/unpaired) as input. See MLFeatures(loop) for details.
	 */
	Array3D<energy_t> MakeStackingTable(const std::vector<StackingFeature> &features) const {
		using namespace std;
		unsigned Nunpaired = 0, Nbranches = 0;
		for (const auto &ft : features) {
			if (ft.i == ft.j)
				++Nunpaired;
			else
				++Nbranches;
		}
		int maxA = Nunpaired + Nbranches, maxB = Nbranches;
		int NF = static_cast<int>(features.size());
		// +3 to maxA and maxB to allow a closing branch to be considered later.
		// Another +3 to avoid needed to code extra cases.
		Array3D<energy_t> dp_table(static_cast<unsigned>(maxA) + 7,
								   static_cast<unsigned>(maxB) + 7,
								   static_cast<unsigned>(NF) + 1,
								   em.MaxMFE());
		// Base cases are an empty suffix with >=0 as and bs.
		// Up to +3 to allow closing branches to be considered.
		for (int a = 0; a <= maxA + 3; ++a) { // a length segments left
			for (int b = 0; b <= maxB + 3; ++b) { // b length segments left
				dp_table[a][b][NF] = em.MLInit(a, b); // Account for a and b.
			}
		}
		for (int i = NF - 1; i >= 0; --i) { // Suffix index.
			for (int a = 0; a <= maxA + 3; ++a) { // a length segments left
				for (int b = 0; b <= maxB + 3; ++b) { // b length segments left=
					int fts_left = NF - i;
					energy_t best = em.MaxMFE();
					// Feature i is not in an interaction.
					if (IsBranch(features[i])) {
						best = dp_table[a + 1][b + 1][i + 1];
					} else if (!IsBranch(features[i])) {
						best = dp_table[a + 1][b][i + 1];
					}
					// Feature i is in some sort of interaction.
					if (fts_left >= 2) {
						// Dangling ends
						// (_).
						if (IsBranch(features[i]) && !IsBranch(features[i + 1])) {
							best =
								min(best, dp_table[a + 2][b + 1][i + 2] + em.ThreeDangle(features[i].i, features[i].j));
						}
						// .(_)
						if (!IsBranch(features[i]) && IsBranch(features[i + 1])) {
							best = min(best,
									   dp_table[a + 2][b + 1][i + 2]
										   + em.FiveDangle(features[i + 1].i, features[i + 1].j));
						}
						// Flush coax (_)(_)
						if (IsBranch(features[i]) && IsBranch(features[i + 1])) {
							best = min(best, dp_table[a + 2][b][i + 2]
								+ em.FlushCoax(features[i].i, features[i].j, features[i + 1].i, features[i + 1].j));
						}
					}
					if (fts_left >= 3) {
						// Terminal mismatch
						// .(_).
						if (!IsBranch(features[i]) && IsBranch(features[i + 1]) && !IsBranch(features[i + 2])) {
							best = min(best,
									   dp_table[a + 3][b + 1][i + 3]
										   + em.Mismatch(features[i + 1].i, features[i + 1].j));
						}
					}
					if (fts_left >= 4) {
						// Mismatch coax
						// (_).(_).
						if (IsBranch(features[i]) && !IsBranch(features[i + 1]) && IsBranch(features[i + 2])
							&& !IsBranch(features[i + 3])) {
							best = min(best, dp_table[a + 2][b][i + 4]
								+ em.MismatchCoax(features[i + 2].i, features[i + 2].j, features[i].i, features[i].j));
						}
						// .(_).(_)
						if (!IsBranch(features[i]) && IsBranch(features[i + 1]) && !IsBranch(features[i + 2])
							&& IsBranch(features[i + 3])) {
							best = min(best, dp_table[a + 2][b][i + 4]
								+ em.MismatchCoax(features[i + 1].i, features[i + 1].j, features[i + 3].i,
												  features[i + 3].j));
						}
					}
					dp_table[a][b][i] = best;
				}
			}
		}
		return dp_table;
	}

public:


	// TODO Override trace functions

	/**
	 * Calculates the FE of the optimal configuration of dangles for a multiloop surface.
	 * Also considers coaxial stacking and a/b length  segments.
	 */
	virtual std::tuple<energy_t, energy_t> OptimalMLConfig(const Surface &surf) const override {
		using namespace std;

		LoopRegion loop(surf);
		auto fts = MLStackingFeatures(loop);
		auto dp_table = MakeStackingTable(fts);
		// Closing branch does nothing.
		energy_t best = dp_table[1][1][0];
		if (!IsBranch(fts[0])) {
			// (._)
			best = min(best, dp_table[2][1][1] + em.ClosingThreeDangle(loop.i, loop.j));
			// (.(_)._)
			if (IsBranch(fts[1]) && !IsBranch(fts[2])) {
				best = min(best, dp_table[2][0][3] + em.MismatchCoax(fts[1].i, fts[1].j, loop.i, loop.j));
			}
		} else { // ((_)_)
			best = min(best, dp_table[2][0][1] + em.FlushCoax(loop.i, loop.j, fts[0].i, fts[0].j));
		}
		// Pop the back.
		auto back = fts.back();
		fts.pop_back();
		auto dptmp = MakeStackingTable(fts);
		if (!IsBranch(back)) {
			// (_.)
			best = min(best, dptmp[2][1][0] + em.ClosingFiveDangle(loop.i, loop.j));
			if (!IsBranch(fts[0])) {
				// (._.)
				best = min(best, dptmp[3][1][1] + em.ClosingMismatch(loop.i, loop.j));
				if (IsBranch(fts.back())) {
					auto back2 = fts.back();
					fts.pop_back();
					auto dptmp2 = MakeStackingTable(fts);
					// (._(_).)
					best = min(best, dptmp2[2][0][1] + em.MismatchCoax(loop.i, loop.j, back2.i, back2.j));
					fts.push_back(back2);
				}
				if (IsBranch(fts[1])) {
					// (.(_)_.)
					best = min(best, dptmp[2][0][2] + em.MismatchCoax(loop.i, loop.j, fts[1].i, fts[1].j));
				}
			}
			if (IsBranch(fts.back()) && !IsBranch(fts[fts.size() - 2])) {
				auto back2 = fts.back();
				fts.pop_back();
				auto back3 = fts.back();
				fts.pop_back();
				auto dptmp3 = MakeStackingTable(fts);
				// (_.(_).)
				best = min(best, dptmp3[2][0][0] + em.MismatchCoax(back2.i, back2.j, loop.i, loop.j));
				fts.push_back(back3);
				fts.push_back(back2);
			}

		} else { // (_(_))
			best = min(best, dptmp[2][0][0] + em.FlushCoax(loop.i, loop.j, back.i, back.j));
		}
		fts.push_back(back);
		// TODO Mkae it return stacking and closure instead of both at once
		return make_tuple(0, best);
	}


	AalbertsScorer(AalbertsModel _em)
		: NNScorer(_em) {}

	AalbertsScorer(AalbertsModel _em, bool _stacking)
		: NNScorer(_em, _stacking) {}
};
}

#endif //RNARK_AALBERTS_SCORER_HPP

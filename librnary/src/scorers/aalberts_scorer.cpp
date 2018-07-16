//
// Created by max on 7/24/17.
//

#include "scorers/aalberts_scorer.hpp"

using namespace librnary;
using namespace std;

Array3D<energy_t> librnary::AalbertsScorer::MakeStackingTable(const vector<StackingFeature> &features) const {
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

vector<Stacking> librnary::AalbertsScorer::TraceStackingTable(const std::vector<StackingFeature> &features,
															  const Array3D<energy_t> &dp_table,
															  int a,
															  int b,
															  int i) const {
	unsigned Nunpaired = 0, Nbranches = 0;
	for (const auto &ft : features) {
		if (ft.i == ft.j)
			++Nunpaired;
		else
			++Nbranches;
	}
	int NF = static_cast<int>(features.size());

	vector<Stacking> stacks;

	while (i < NF) {
		int nexta = -1, nextb = -1, nexti = -1, fts_left = NF - i;
		librnary::energy_t best = em.MaxMFE(), v;
		Stacking stack(NONE, -1, -1);
		// Feature i is not in an interaction.
		if (IsBranch(features[i])) {
			tie(nexta, nextb, nexti) = make_tuple(a + 1, b + 1, i + 1);
			best = dp_table[a + 1][b + 1][i + 1];
		} else if (!IsBranch(features[i])) {
			tie(nexta, nextb, nexti) = make_tuple(a + 1, b, i + 1);
			best = dp_table[a + 1][b][i + 1];
		}
		// Feature i is in some sort of interaction.
		if (fts_left >= 2) {
			// Dangling ends
			// (_).
			if (IsBranch(features[i]) && !IsBranch(features[i + 1])) {
				v = dp_table[a + 2][b + 1][i + 2] + em.ThreeDangle(features[i].i, features[i].j);
				if (v < best) {
					best = v;
					tie(nexta, nextb, nexti) = make_tuple(a + 2, b + 1, i + 2);
					stack = Stacking(DANLGE3, features[i].i, features[i].j);
				}
			}
			// .(_)
			if (!IsBranch(features[i]) && IsBranch(features[i + 1])) {
				v = dp_table[a + 2][b + 1][i + 2] + em.FiveDangle(features[i + 1].i, features[i + 1].j);
				if (v < best) {
					best = v;
					tie(nexta, nextb, nexti) = make_tuple(a + 2, b + 1, i + 2);
					stack = Stacking(DANLGE5, features[i + 1].i, features[i + 1].j);
				}
			}
			// Flush coax (_)(_)
			if (IsBranch(features[i]) && IsBranch(features[i + 1])) {
				v = dp_table[a + 2][b][i + 2]
					+ em.FlushCoax(features[i].i, features[i].j, features[i + 1].i, features[i + 1].j);
				if (v < best) {
					best = v;
					tie(nexta, nextb, nexti) = make_tuple(a + 2, b, i + 2);
					stack = Stacking(FLUSHCX, features[i].i, features[i].j, features[i + 1].i, features[i + 1].j);
				}
			}
		}
		if (fts_left >= 3) {
			// Terminal mismatch
			// .(_).
			if (!IsBranch(features[i]) && IsBranch(features[i + 1]) && !IsBranch(features[i + 2])) {
				v = dp_table[a + 3][b + 1][i + 3] + em.Mismatch(features[i + 1].i, features[i + 1].j);
				if (v < best) {
					best = v;
					tie(nexta, nextb, nexti) = make_tuple(a + 3, b + 1, i + 3);
					stack = Stacking(MISMATCH, features[i + 1].i, features[i + 1].j);
				}
			}
		}
		if (fts_left >= 4) {
			// Mismatch coax
			// (_).(_).
			if (IsBranch(features[i]) && !IsBranch(features[i + 1]) && IsBranch(features[i + 2])
				&& !IsBranch(features[i + 3])) {
				v = dp_table[a + 2][b][i + 4]
					+ em.MismatchCoax(features[i + 2].i, features[i + 2].j, features[i].i, features[i].j);
				if (v < best) {
					best = v;
					tie(nexta, nextb, nexti) = make_tuple(a + 2, b, i + 4);
					stack = Stacking(MMCX, features[i + 2].i, features[i + 2].j, features[i].i, features[i].j);
				}
			}
			// .(_).(_)
			if (!IsBranch(features[i]) && IsBranch(features[i + 1]) && !IsBranch(features[i + 2])
				&& IsBranch(features[i + 3])) {
				v = dp_table[a + 2][b][i + 4]
					+ em.MismatchCoax(features[i + 1].i, features[i + 1].j, features[i + 3].i, features[i + 3].j);
				if (v < best) {
					best = v;
					tie(nexta, nextb, nexti) = make_tuple(a + 2, b, i + 4);
					stack = Stacking(MMCX, features[i + 1].i, features[i + 1].j, features[i + 3].i, features[i + 3].j);
				}
			}
		}
		a = nexta;
		b = nextb;
		i = nexti;
		if (stack.t != NONE)
			stacks.push_back(stack);
	}
	return stacks;

}

std::tuple<energy_t, energy_t> librnary::AalbertsScorer::OptimalMLConfig(const Surface &surf) const {
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
	// TODO Make it return stacking and closure instead of both at once
	return make_tuple(0, best);
}

std::tuple<std::vector<Stacking>, ClosureFeaturesMap>
librnary::AalbertsScorer::TraceMLConfig(const Surface &surf) const {
	LoopRegion loop(surf);
	auto fts = MLStackingFeatures(loop);
	auto dp_table = MakeStackingTable(fts);
	// Closing branch does nothing.
	vector<Stacking> stacks = TraceStackingTable(fts, dp_table,1, 1, 0);
	energy_t best = dp_table[1][1][0], v;
	if (!IsBranch(fts[0])) {
		// (._)
		v = dp_table[2][1][1] + em.ClosingThreeDangle(loop.i, loop.j);
		if (v < best) {
			best = v;
			stacks = TraceStackingTable(fts, dp_table, 2, 1, 1);
			stacks.emplace_back(DANLGE3, loop.i, loop.j);
		}
		// (.(_)._)
		if (IsBranch(fts[1]) && !IsBranch(fts[2])) {
			v = dp_table[2][0][3] + em.MismatchCoax(fts[1].i, fts[1].j, loop.i, loop.j);
			if (v < best) {
				best = v;
				stacks = TraceStackingTable(fts, dp_table,2, 0, 3);
				stacks.emplace_back(MMCX, fts[1].i, fts[1].j, loop.i, loop.j);
			}
		}
	} else { // ((_)_)
		v = dp_table[2][0][1] + em.FlushCoax(loop.i, loop.j, fts[0].i, fts[0].j);
		if (v < best) {
			best = v;
			stacks = TraceStackingTable(fts, dp_table,2, 0, 1);
			stacks.emplace_back(FLUSHCX, loop.i, loop.j, fts[0].i, fts[0].j);
		}
	}
	// Pop the back.
	auto back = fts.back();
	fts.pop_back();
	auto dptmp = MakeStackingTable(fts);
	if (!IsBranch(back)) {
		// (_.)
		v = dptmp[2][1][0] + em.ClosingFiveDangle(loop.i, loop.j);
		if (v < best) {
			best = v;
			stacks = TraceStackingTable(fts, dptmp, 2, 1, 0);
			stacks.emplace_back(DANLGE5, loop.i, loop.j);
		}
		if (!IsBranch(fts[0])) {
			// (._.)
			v = dptmp[3][1][1] + em.ClosingMismatch(loop.i, loop.j);
			if (v < best) {
				best = v;
				stacks = TraceStackingTable(fts, dptmp, 3, 1, 1);
				stacks.emplace_back(MISMATCH, loop.i, loop.j);
			}
			if (IsBranch(fts.back())) {
				auto back2 = fts.back();
				fts.pop_back();
				auto dptmp2 = MakeStackingTable(fts);
				// (._(_).)
				v = dptmp2[2][0][1] + em.MismatchCoax(loop.i, loop.j, back2.i, back2.j);
				if (v < best) {
					best = v;
					stacks = TraceStackingTable(fts, dptmp2, 2,0,1);
					stacks.emplace_back(MMCX, loop.i, loop.j, back2.i, back2.j);
				}
				fts.push_back(back2);
			}
			if (IsBranch(fts[1])) {
				// (.(_)_.)
				v = dptmp[2][0][2] + em.MismatchCoax(loop.i, loop.j, fts[1].i, fts[1].j);
				if (v < best) {
					best = v;
					stacks = TraceStackingTable(fts, dptmp, 2,0,2);
					stacks.emplace_back(MMCX, loop.i, loop.j, fts[1].i, fts[1].j);
				}
			}
		}
		if (IsBranch(fts.back()) && !IsBranch(fts[fts.size() - 2])) {
			auto back2 = fts.back();
			fts.pop_back();
			auto back3 = fts.back();
			fts.pop_back();
			auto dptmp3 = MakeStackingTable(fts);
			// (_.(_).)
			v = dptmp3[2][0][0] + em.MismatchCoax(back2.i, back2.j, loop.i, loop.j);
			if (v < best) {
				best = v;
				stacks = TraceStackingTable(fts, dptmp3, 2,0,0);
				stacks.emplace_back(MMCX, back2.i, back2.j, loop.i, loop.j);
			}
			fts.push_back(back3);
			fts.push_back(back2);
		}

	} else { // (_(_))
		v = dptmp[2][0][0] + em.FlushCoax(loop.i, loop.j, back.i, back.j);
		if (v < best) {
			best = v;
			stacks = TraceStackingTable(fts, dptmp, 2,0,0);
			stacks.emplace_back(FLUSHCX, loop.i, loop.j, back.i, back.j);
		}
	}
	fts.push_back(back);

	int lengthA = 0, lengthB = 0;
	tie(lengthA, lengthB) = ExtractLengthALengthB(stacks, loop);


	return std::tuple<std::vector<Stacking>, ClosureFeaturesMap>(stacks, {{"N", lengthA}, {"M", lengthB}});
}
//
// Created by max on 6/14/16.
//

#ifndef RNARK_NN_SCORER_HPP
#define RNARK_NN_SCORER_HPP

#include <string>
#include <sstream>
#include <vector>
#include <stack>
#include <tuple>
#include <algorithm>
#include <map>

#include "ss_tree.hpp"
#include "multi_array.hpp"
#include "multi_loop.hpp"
#include "energy.hpp"

namespace librnary {

typedef std::map<std::string, energy_t> ClosureFeaturesMap;

namespace {
std::string ToString(const Surface &surf) {
	return "(" + std::to_string(surf.PairI()) + ", " + std::to_string(surf.PairJ()) + ")";
}
}

/**
 * A simple tree data structure that represents that energy score infomation for a secondary structure surface.
 */
struct SurfaceScore {
	Surface surface;
	librnary::energy_t loop_score = 0, recursive_score = 0;
	librnary::energy_t AUGUclosure = 0, stacking = 0, ml_closure;
	ClosureFeaturesMap ml_closure_features;
	std::vector<Stacking> stacks;
	std::vector<std::shared_ptr<SurfaceScore>> subsurfaces;
	SurfaceScore(const Surface &surf) : surface(surf), loop_score(0), recursive_score(0),
										AUGUclosure(0), stacking(0), ml_closure(0) {}
	virtual std::string Describe(char indent_char, unsigned indent_level) const {
		std::stringstream ss;
		std::string indent(indent_level, indent_char);
		ss << indent;
		if (surface.IsExternalLoop()) {
			ss << "External-loop: " << loop_score << "/" << recursive_score << std::endl;
			ss << indent << "AU/GU penalties: " << AUGUclosure << std::endl;
			ss << indent << "Stacking: " << stacking << std::endl;
			ss << indent;
			for (const auto &s : stacks) {
				ss << s.ToString() << " ";
			}
			ss << std::endl;
			for (const auto &subsurf : subsurfaces) {
				ss << subsurf->Describe(indent_char, indent_level + 1);
			}

		} else if (surface.NumChildren() == 0) {
			ss << "One-loop closed by " << ToString(surface)
			   << ": " << loop_score << "/" << recursive_score << std::endl;
		} else if (surface.NumChildren() == 1) {
			ss << "Two-loop closed by " << ToString(surface) << " and " << ToString(surface.Child(0))
			   << ": " << loop_score << "/" << recursive_score << std::endl;
			ss << subsurfaces[0]->Describe(indent_char, indent_level + 1);
		} else {
			ss << "Multi-loop closed by " << ToString(surface) << ": "
			   << loop_score << "/" << recursive_score << std::endl;
			ss << indent << "AU/GU penalties: " << AUGUclosure << std::endl;
			ss << indent << "Multi-loop closure: " << ml_closure << std::endl;
			ss << indent << "Multi-loop closure featues: ";
			for (const auto &kv : ml_closure_features) {
				ss << kv.first << "=" << kv.second << " ";
			}
			ss << std::endl;
			ss << indent << "Stacking: " << stacking << std::endl;
			ss << indent;
			for (const auto &s : stacks) {
				ss << s.ToString() << " ";
			}
			ss << std::endl;
			for (const auto &subsurf : subsurfaces) {
				ss << subsurf->Describe(indent_char, indent_level + 1);
			}
		}
		return ss.str();
	}

};

template<typename EnergyModel>
/**
 * Given a sequence and a structure, computes the free energy score of the structure.
 * Also implements related functionality.
 * Essentially an "efn" class.
 */
class NNScorer {
private:
	// Trace back through the stacking table to find the optimal stacking.
	// Returns a list of stack types and bonding pairs.
	std::vector<Stacking> TraceStacking(const LoopRegion &loop,
										const Array3D<energy_t> &table, int start_push,
										int start_pull, int start_s) const {
		using namespace std;
		int NS = static_cast<int>(loop.enclosed.size());
		stack<tuple<int, int, int>> trace_stack;
		trace_stack.push(make_tuple(start_push, start_pull, start_s));
		vector<Stacking> output;
		while (!trace_stack.empty()) {
			tuple<int, int, int> state = trace_stack.top();
			trace_stack.pop();
			int push = get<0>(state), pull = get<1>(state), s = get<2>(state);
			// Base cases.
			if (s == NS)
				continue;
			int i = loop.i + push, j = loop.j - pull;
			if (s != 0)
				i = loop.enclosed[s - 1].j + push;
			if (s != NS - 1)
				j = loop.enclosed[s + 1].i;
			auto ss = loop.enclosed[s];
			if (ss.i <= i)
				continue;
			// Recurse normally.
			energy_t rec = table[0][pull][s + 1];
			// Recurse with extra padding on the right.
			energy_t recr = table[1][pull][s + 1];
			// Find min(no dangle, 3' dangle, 5' dangle, terminal em.Mismatch).
			// rec and recr captures recursively dealing with 3' branches.
			energy_t best = rec;
			Stacking decomp(NONE, ss.i, ss.j);
			tuple<int, int, int> rec_n(0, pull, s + 1);
			tuple<int, int, int> recr_n(1, pull, s + 1);
			auto next = rec_n;
			if (ss.j + 1 < j) {
				energy_t score = recr + em.ThreeDangle(ss.i, ss.j);
				if (score < best) {
					best = score;
					decomp = Stacking(DANLGE3, ss.i, ss.j);
					next = recr_n;
				}
			}
			if (i + 1 < ss.i) {
				energy_t score = rec + em.FiveDangle(ss.i, ss.j);
				if (score < best) {
					best = score;
					decomp = Stacking(DANLGE5, ss.i, ss.j);
					next = rec_n;
				}
				if (ss.j + 1 < j) {
					score = recr + em.Mismatch(ss.i, ss.j);
					if (score < best) {
						best = score;
						decomp = Stacking(MISMATCH, ss.i, ss.j);
						next = recr_n;
					}
				}
			}

			// Coaxial stacks.
			if (s + 1 < NS) { // If there is a subsequent branch.
				rec = table[0][pull][s + 2];
				recr = table[1][pull][s + 2];
				rec_n = tuple<int, int, int>(0, pull, s + 2);
				recr_n = tuple<int, int, int>(1, pull, s + 2);
				auto n = loop.enclosed[s + 1];
				if (n.i == ss.j + 1) { // Next branch is flush.
					energy_t score = rec + em.FlushCoax(ss.i, ss.j, n.i, n.j);
					if (score < best) {
						best = score;
						decomp = Stacking(FLUSHCX, ss.i, ss.j, n.i, n.j);
						next = rec_n;
					}
				}
				if (n.i == ss.j + 2) { // There is an unpaired between this and next branch.
					if (i + 1 < ss.i) { // em.Mismatch for this branch.
						energy_t score = rec + em.MismatchCoax(ss.i, ss.j, n.i, n.j);
						if (score < best) {
							best = score;
							decomp = Stacking(MMCX, ss.i, ss.j, n.i, n.j);
							next = rec_n;
						}
					}
					if (n.j + 1 < loop.j - pull) {
						energy_t score = recr + em.MismatchCoax(n.i, n.j, ss.i, ss.j);
						if (score < best) {
							// No need to assign best as it won't be used. I've left the line in here as a reminder.
							// Just in cast this code gets copied and pasted.
							// best = score;
							decomp = Stacking(MMCX, n.i, n.j, ss.i, ss.j);
							next = recr_n;
						}
					}
				}
			}

			trace_stack.push(next);
			output.push_back(decomp);
		}
		return output;
	}

	// Returns a DP table containing the optimal configurations of branches.
	// dp_table[i][j][s] is the MFE configuration of stacking from subsurface s to NS-1.
	// The immediate left subsurface used i right dangles, and the closing surface used j on the right perimeter.
	// Optimal "stacking" is defined here http://rna.urmc.rochester.edu/NNDB/turner04/exterior.html
	Array3D<energy_t> MakeStackingTable(const LoopRegion &loop) const {
		using namespace std;
		int NS = static_cast<int>(loop.enclosed.size());
		Array3D<energy_t> dp_table(2, 2, static_cast<unsigned>(NS + 1), em.MaxMFE());
//		vector<vector<vector<int>>> dp_table(2, vector<vector<int>>(2, vector<int>(static_cast<unsigned>(NS + 1),
//																				   em.MaxMFE())));
		dp_table[0][0][NS] = 0;
		if (loop.enclosed[NS - 1].j + 1 < loop.j) {
			dp_table[1][0][NS] = 0;
			dp_table[0][1][NS] = 0;
		}
		if (loop.enclosed[NS - 1].j + 1 < loop.j - 1)
			dp_table[1][1][NS] = 0;
		for (int s = NS - 1; s >= 0; --s) {
			for (int push = 0; push < 2; ++push) {
				for (int pull = 0; pull < 2; ++pull) {
					// i and j are the left and right limits for subsurface s.
					int i = loop.i + push, j = loop.j - pull;
					if (s != 0)
						i = loop.enclosed[s - 1].j + push;
					if (s != NS - 1)
						j = loop.enclosed[s + 1].i;
					auto ss = loop.enclosed[s];
					if (ss.i <= i)
						continue;
					// Recurse normally.
					energy_t rec = dp_table[0][pull][s + 1];
					// Recurse with extra padding on the right.
					energy_t recr = dp_table[1][pull][s + 1];
					// Find min(no dangle, 3' dangle, 5' dangle, terminal em.Mismatch).
					// rec and recr captures recursively dealing with 3' branches.
					energy_t best = rec;
					if (ss.j + 1 < j)
						best = min(best, recr + em.ThreeDangle(ss.i, ss.j));
					if (i + 1 < ss.i) {
						best = min(best, rec + em.FiveDangle(ss.i, ss.j));
						if (ss.j + 1 < j)
							best = min(best, recr + em.Mismatch(ss.i, ss.j));
					}

					// Coaxial stacks.
					if (s + 1 < NS) { // If there is a subsequent branch.
						rec = dp_table[0][pull][s + 2];
						recr = dp_table[1][pull][s + 2];
						auto next = loop.enclosed[s + 1];
						if (next.i == ss.j + 1) // Next branch is flush.
							best = min(best, rec + em.FlushCoax(ss.i, ss.j, next.i, next.j));
						if (next.i == ss.j + 2) { // There is an unpaired between this and next branch.
							if (i + 1 < ss.i) // em.Mismatch for this branch.
								best = min(best, rec + em.MismatchCoax(ss.i, ss.j, next.i, next.j));
							if (next.j + 1 < loop.j - pull)
								best = min(best, recr + em.MismatchCoax(next.i, next.j, ss.i, ss.j));
						}
					}

					dp_table[push][pull][s] = best;
				}
			}
		}
		return dp_table;
	}
protected:
	EnergyModel em;

	bool stacking = true;

public:

	std::vector<Stacking> TraceExternalStacking(const Surface &surf) const {
		const int NS = surf.NumChildren();
		if (NS > 0) {
			auto dp_table = MakeStackingTable(surf);
			return TraceStacking(LoopRegion(surf), dp_table, 0, 0, 0);
		} else {
			return {};
		}
	}

	/// Calculates the FE of the optimal configuration of dangles for a multiloop surface
	/// Assuming no internal dangles are allowed. This is for external multiloops.
	energy_t OptimalExternalStacking(const Surface &surf) const {
		const int NS = surf.NumChildren();
		if (NS > 0) {
			auto dp_table = MakeStackingTable(surf);
			return dp_table[0][0][0];
		}
		return 0;
	}

	/// Score an exterior super-surface. Is recursive.
	energy_t ScoreExterior(const librnary::Surface &super) const {
		assert(super.IsExternalLoop());
		energy_t sum = 0;
		if (stacking)
			sum += OptimalExternalStacking(super);
		for (const auto &ss : super.Children()) {
			sum += ScoreInternal(ss) + em.Branch(ss.PairI(), ss.PairJ());
		}
		return sum;
	}

	SurfaceScore TraceExterior(const librnary::Surface &super) const {
		SurfaceScore surfscore(super);
		if (stacking) {
			surfscore.stacking = OptimalExternalStacking(super);
			surfscore.stacks = TraceExternalStacking(super);
		}
		for (const auto &ss : super.Children()) {
			surfscore.AUGUclosure += em.Branch(ss.PairI(), ss.PairJ());
			surfscore.recursive_score += ScoreInternal(ss);
			surfscore.subsurfaces.push_back(std::unique_ptr<SurfaceScore>(new SurfaceScore(TraceInternal(ss))));
		}
		surfscore.loop_score = surfscore.stacking + surfscore.AUGUclosure;
		surfscore.recursive_score += surfscore.loop_score;
		assert(surfscore.recursive_score == ScoreExterior(super));
		return surfscore;
	}

	/// Scores an internal surface. Is recursive.
	virtual energy_t ScoreInternal(const librnary::Surface &surf) const {
		assert(surf.PairI() >= 0 && static_cast<int>(em.RNA().size()) > surf.PairJ());
		if (surf.NumChildren() == 0)
			return em.OneLoop(surf.PairI(), surf.PairJ());
		if (surf.NumChildren() == 1)
			return em.TwoLoop(surf.PairI(), surf.Child(0).PairI(), surf.Child(0).PairJ(), surf.PairJ()) +
				ScoreInternal(surf.Child(0));
		// The remaining code considers internal multiloops.
		// Start with closure and closing branch penalty.
		auto mlConfig = OptimalMLConfig(surf);
		energy_t sum = std::get<1>(mlConfig) + em.Branch(surf.PairI(), surf.PairJ());
		if (stacking)
			sum += std::get<0>(OptimalMLConfig(surf));
		for (const auto &ss : surf.Children())
			sum += ScoreInternal(ss) + em.Branch(ss.PairI(), ss.PairJ());
		return sum;
	}

	virtual SurfaceScore TraceInternal(const librnary::Surface &surf) const {
		SurfaceScore sscore(surf);
		if (surf.NumChildren() == 0) {
			sscore.loop_score = em.OneLoop(surf.PairI(), surf.PairJ());
			sscore.recursive_score = sscore.loop_score;
		} else if (surf.NumChildren() == 1) {
			sscore.loop_score = em.TwoLoop(surf.PairI(), surf.Child(0).PairI(), surf.Child(0).PairJ(), surf.PairJ());
			sscore.recursive_score = sscore.loop_score + ScoreInternal(surf.Child(0));
			sscore.subsurfaces.push_back(std::unique_ptr<SurfaceScore>(new SurfaceScore(TraceInternal(surf.Child(0)))));
		} else {
			// The remaining code considers internal multiloops.
			// Start with closure and closing branch penalty.
			auto ml_conf_score = OptimalMLConfig(surf);
			sscore.ml_closure = std::get<1>(ml_conf_score);
			sscore.loop_score = sscore.ml_closure;
			if (stacking) {
				sscore.stacking = std::get<0>(ml_conf_score);
				sscore.loop_score += sscore.stacking;
			}
			std::tuple<std::vector<Stacking>, ClosureFeaturesMap> ml_conf = TraceMLConfig(surf);
			if (stacking)
				sscore.stacks = std::get<0>(ml_conf);
			sscore.ml_closure_features = std::get<1>(ml_conf);
			sscore.AUGUclosure += em.Branch(surf.PairI(), surf.PairJ());
			for (const auto &ss : surf.Children()) {
				sscore.recursive_score += ScoreInternal(ss);
				sscore.AUGUclosure += em.Branch(ss.PairI(), ss.PairJ());
				sscore.subsurfaces.push_back(std::unique_ptr<SurfaceScore>(new SurfaceScore(TraceInternal(ss))));
			}
			sscore.loop_score += sscore.AUGUclosure;
			sscore.recursive_score += sscore.loop_score;
		}
		assert(sscore.recursive_score == ScoreInternal(surf));
		return sscore;
	}

	/// Calculates the FE of the optimal configuration of dangles for a multiloop.
	/// Returns a tuple of stacking energy and closure energy.
	virtual energy_t OptimalMLStacking(const LoopRegion &loop) const {
		using namespace std;


		auto dp_table = MakeStackingTable(loop);
		// Closing helix has a dangle or terminal em.Mismatch.
		energy_t best = min(
			dp_table[1][0][0] + em.ClosingThreeDangle(loop.i, loop.j),
			dp_table[0][1][0] + em.ClosingFiveDangle(loop.i, loop.j));
		best = min(best, dp_table[1][1][0] + em.ClosingMismatch(loop.i, loop.j));
		best = min(best, dp_table[0][0][0]); // No dangle/em.Mismatch, just plain helix.
		// Coaxial stacking.
		// Note that in the if statements we can safely assume that loop.enclosed has size >= 2.
		// This is because we are definitely in a multiloop, and multiloops have >= 2 internal branches.
		// ((_)_)
		if (loop.enclosed[0].i == loop.i + 1)
			best = min(best, dp_table[0][0][1]
				+ em.FlushCoax(loop.i, loop.j, loop.enclosed[0].i, loop.enclosed[0].j));

		if (loop.enclosed[0].i == loop.i + 2) {
			// (.(_)._)
			if (loop.enclosed[1].i > loop.enclosed[0].j + 1)
				best = min(best, dp_table[1][0][1]
					+ em.MismatchCoax(loop.enclosed[0].i, loop.enclosed[0].j, loop.i, loop.j));
			// (.(_)_.)
			if (loop.enclosed[0].i == loop.i + 2)
				best = min(best, dp_table[0][1][1] +
					em.MismatchCoax(loop.i, loop.j, loop.enclosed[0].i, loop.enclosed[0].j));
		}

		// The following lines are going to seem crazy.
		// But they are actually an elegant (if I do say so myself) way to avoid the O(N^2) DP.
		// It also avoids coding a right to left version as well as the left to right version.
		// This way we can simply execute the O(N) DP twice.
		auto final_helix = loop.enclosed[loop.enclosed.size() - 1];
		LoopRegion loop_coaxright = loop;
		loop_coaxright.j = final_helix.i;
		loop_coaxright.enclosed.pop_back();
		auto dp_table_Cxright = MakeStackingTable(loop_coaxright);
		// (_(_))
		if (final_helix.j == loop.j - 1)
			best = min(best, dp_table_Cxright[0][0][0]
				+ em.FlushCoax(loop.i, loop.j, final_helix.i, final_helix.j));
		// (_.(_).)
		if (final_helix.j == loop.j - 2) {
			if (loop.enclosed[loop.enclosed.size() - 2].j < final_helix.i - 1)
				best = min(best, dp_table_Cxright[0][1][0]
					+ em.MismatchCoax(final_helix.i, final_helix.j, loop.i, loop.j));
			// (._(_).)
			best = min(best,
					   dp_table_Cxright[1][0][0] + em.MismatchCoax(loop.i, loop.j, final_helix.i, final_helix.j));

		}
		return best;
	}

	/**
	 * Computes the score of optimal ML stacking configuration.
	 * @param surf The multi-loop.
	 * @return A tuple of optimal ML stacking score and ML closure score.
	 */
	virtual std::tuple<energy_t, energy_t> OptimalMLConfig(const Surface &surf) const {
		return std::tuple<energy_t, energy_t>(OptimalMLStacking(surf), em.MLClosure(surf));
	}


	/// Returns a tuple of stacking interation list, and closure features.
	virtual std::tuple<std::vector<Stacking>, ClosureFeaturesMap> TraceMLConfig(const Surface &surf) const {
		using namespace std;

		LoopRegion loop(surf);

		auto dp_table = MakeStackingTable(LoopRegion(surf));
		// Normal.
		vector<Stacking> stacks = TraceStacking(surf, dp_table, 0, 0, 0);
		stacks.emplace_back(NONE, loop.i, loop.j);
		energy_t best = dp_table[0][0][0];
		// 3' dangle.
		energy_t score = dp_table[1][0][0] + em.ClosingThreeDangle(loop.i, loop.j);
		if (score < best) {
			best = score;
			stacks = TraceStacking(surf, dp_table, 1, 0, 0);
			stacks.emplace_back(DANLGE3, loop.i, loop.j);
		}
		// 5' dangle.
		score = dp_table[0][1][0] + em.ClosingFiveDangle(loop.i, loop.j);
		if (score < best) {
			best = score;
			stacks = TraceStacking(surf, dp_table, 0, 1, 0);
			stacks.emplace_back(DANLGE5, loop.i, loop.j);
		}
		// Mismatch
		score = dp_table[1][1][0] + em.ClosingMismatch(loop.i, loop.j);
		if (score < best) {
			best = score;
			stacks = TraceStacking(surf, dp_table, 1, 1, 0);
			stacks.emplace_back(MISMATCH, loop.i, loop.j);
		}
		// Coaxial stacking.
		// Note that in the if statements we can safely assume that loop.enclosed has size >= 2.
		// This is because we are definitely in a multiloop, and multiloops have >= 2 internal branches.
		// ((_)_)
		if (loop.enclosed[0].i == loop.i + 1) {
			score = dp_table[0][0][1] +
				em.FlushCoax(loop.i, loop.j, loop.enclosed[0].i, loop.enclosed[0].j);
			if (score < best) {
				best = score;
				stacks = TraceStacking(surf, dp_table, 0, 0, 1);
				stacks.emplace_back(FLUSHCX, loop.i, loop.j, loop.enclosed[0].i, loop.enclosed[0].j);
			}
		}

		if (loop.enclosed[0].i == loop.i + 2) {
			// (.(_)._)
			if (loop.enclosed[1].i > loop.enclosed[0].j + 1) {
				score = dp_table[1][0][1] +
					em.MismatchCoax(loop.enclosed[0].i, loop.enclosed[0].j, loop.i, loop.j);
				if (score < best) {
					best = score;
					stacks = TraceStacking(surf, dp_table, 1, 0, 1);
					stacks.emplace_back(MMCX, loop.enclosed[0].i, loop.enclosed[0].j, loop.i, loop.j);
				}
			}
			// (.(_)_.)
			if (loop.enclosed[0].i == loop.i + 2) {
				score = dp_table[0][1][1] +
					em.MismatchCoax(loop.i, loop.j, loop.enclosed[0].i, loop.enclosed[0].j);
				if (score < best) {
					best = score;
					stacks = TraceStacking(surf, dp_table, 0, 1, 1);
					stacks.emplace_back(MMCX, loop.i, loop.j, loop.enclosed[0].i, loop.enclosed[0].j);
				}
			}
		}

		// The following lines are going to seem crazy.
		// But they are actually an elegant (if I do say so myself) way to avoid the O(N^2) DP.
		// It also avoids coding a right to left version as well as the left to right version.
		// This way we can simply execute the O(N) DP twice.
		auto final_helix = loop.enclosed[loop.enclosed.size() - 1];
		LoopRegion loop_coaxright = loop;
		loop_coaxright.j = final_helix.i;
		loop_coaxright.enclosed.pop_back();
		auto dp_table_Cxright = MakeStackingTable(loop_coaxright);
		// (_(_))
		if (final_helix.j == loop.j - 1) {
			score = dp_table_Cxright[0][0][0] + em.FlushCoax(loop.i, loop.j, final_helix.i, final_helix.j);
			if (score < best) {
				best = score;
				stacks = TraceStacking(loop_coaxright, dp_table_Cxright, 0, 0, 0);
				stacks.emplace_back(FLUSHCX, loop.i, loop.j, final_helix.i, final_helix.j);
			}
		}
		// (_.(_).)
		if (final_helix.j == loop.j - 2) {
			if (loop.enclosed[loop.enclosed.size() - 2].j < final_helix.i - 1) {
				score = dp_table_Cxright[0][1][0] + em.MismatchCoax(final_helix.i, final_helix.j, loop.i, loop.j);
				if (score < best) {
					best = score;
					stacks = TraceStacking(loop_coaxright, dp_table_Cxright, 0, 1, 0);
					stacks.emplace_back(MMCX, final_helix.i, final_helix.j, loop.i, loop.j);
				}
			}
			// (._(_).)
			score = dp_table_Cxright[1][0][0] + em.MismatchCoax(loop.i, loop.j, final_helix.i, final_helix.j);
			if (score < best) {
				best = score;
				stacks = TraceStacking(loop_coaxright, dp_table_Cxright, 1, 0, 0);
				stacks.emplace_back(MMCX, loop.i, loop.j, final_helix.i, final_helix.j);
			}
		}
		assert(best == std::get<0>(OptimalMLConfig(surf)));
		ClosureFeaturesMap cfm = {{"unpaired", surf.Unpaired()}, {"branches", surf.NumChildren() + 1}};
		return make_tuple(stacks, cfm);
	}

	void SetRNA(const librnary::PrimeStructure &_rna) {
		em.SetRNA(_rna);
	}

	int MaxMFE() const {
		return em.MaxMFE();
	}

	EnergyModel GetEnergyModel() const {
		return this->em;
	}

	explicit NNScorer(EnergyModel _em)
		: em(_em) {}

	NNScorer(EnergyModel _em, bool _stacking)
		: em(_em), stacking(_stacking) {}

	bool GetStacking() const {
		return stacking;
	}

	void SetStacking(bool v) {
		NNScorer::stacking = v;
	}

};
}

#endif //RNARK_NN_SCORER_HPP

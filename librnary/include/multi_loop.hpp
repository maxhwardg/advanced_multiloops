//
// Created by max on 8/11/16.
//

#ifndef RNARK_MULTI_LOOP_HPP
#define RNARK_MULTI_LOOP_HPP

#include "ss_tree.hpp"

namespace librnary {


enum StackType {
	NONE = 0, DANLGE5, DANLGE3, MISMATCH, FLUSHCX, MMCX
};

/// A stacking interaction.
struct Stacking {
	StackType t;
	int i, j, k, l;
	Stacking(StackType _t, int _i, int _j)
		: t(_t), i(_i), j(_j), k(0), l(0) {}

	Stacking(StackType _t, int _i, int _j, int _k, int _l)
		: t(_t), i(_i), j(_j), k(_k), l(_l) {};

	std::string ToString() const;
};

struct LoopRegion {
	int i, j;
	std::vector<librnary::BondPair> enclosed;
	LoopRegion(const Surface &surf)
		: i(surf.PairI()), j(surf.PairJ()) {
		enclosed.assign(static_cast<unsigned>(surf.NumChildren()), librnary::BondPair());
		for (int idx = 0; idx < surf.NumChildren(); ++idx) {
			enclosed[idx].i = surf.Child(idx).PairI();
			enclosed[idx].j = surf.Child(idx).PairJ();
		}
	}
	LoopRegion(int _i, int _j) : i(_i), j(_j) {}
};

/**
 * A multi-loop (or external-loop) feature that could be involved in stacking.
 * Either a branch or unpaired.
 */
struct StackingFeature {
	int i, j;
	bool IsBranch() const;
	bool IsUnpaired() const;
};

bool IsBranch(const StackingFeature &f);

bool IsUnpaired(const StackingFeature &f);

/**
 * Returns a list of ML features given a multi-loop.
 * This consists of all branches and unpaired nucletoides.
 * Does not consider the closing loop as a branch.
 * A Feature with i==j represents an unpaired/free nucleotide feature.
 */
std::vector<StackingFeature> MLStackingFeatures(const LoopRegion &loop);

/**
 * @param surf The multi-loop surface.
 * @return The sum of asymmetry values for every branch (including closing).
 */
int SumAsymmetry(const Surface &surf);

/**
 * @param num_nucleotides The number of nucleotides incident to the loop. Might be paired or unpaired.
 * @param min_hairpin_unpaired The minimum size of unpaired nucleotides allowed in a hairpin loop.
 * @param branches The number of branches in the multi-loop. If this could be anything, set to 3 as this is the "worst".
 * @return The maximum SumAsymmetry of any multi-loop configuration.
 */
int SumAsymmetryUpperBound(unsigned num_nucleotides, unsigned min_hairpin_unpaired, unsigned branches = 3);


/**
 * Extracts a list of multi-loops from an RNA structure.
 * @param mls The list of multi-loops. This list will have the multi-loops appended to it.
 * @param surf The structure.
 */
void ExtractMultiLoops(std::vector<librnary::LoopRegion> &mls, const librnary::Surface &surf);

/**
 * @param loop A multi-loop.
 * @return The number of branches in the multi-loop.
 */
int ExtractBranches(const librnary::LoopRegion &loop);

/**
 * Extracts unpaired segment lengths from a multi-loop. See asymmetry folders to understand what an unpaired segment is.
 * @param loop The multi-loop.
 * @return A list of the unpaired gap/segment lengths in the multi-loop.
 */
std::vector<int> ExtractSegLengths(const librnary::LoopRegion &loop);

/**
 * Extracts the number of unpaired in a multi-loop.
 * @param loop The multi-loop.
 * @return The number of unpaired nts.
 */
int ExtractUnpaired(const librnary::LoopRegion &loop);

/**
 * Extracts the sum of asymmetries from a multi-loop.
 * @param loop The multi-loop.
 * @return The sum of branch "asymmetries" according the the Matthews 2002 definition.
 */
int ExtractSumAsymmetry(const librnary::LoopRegion &loop);

}

#endif //RNARK_MULTI_LOOP_HPP

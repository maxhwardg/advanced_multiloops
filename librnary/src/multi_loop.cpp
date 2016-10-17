//
// Created by max on 8/11/16.
//

#include "multi_loop.hpp"

#include <sstream>

using namespace std;

int librnary::SumAsymmetry(const librnary::Surface &surf) {
	assert(surf.NumChildren() >= 2);
	int sum_asymmetry = 0;
	int first_up = surf.Child(0).PairI() - surf.PairI() - 1;
	int prev_up = first_up;
	for (int i = 0; i < surf.NumChildren() - 1; ++i) {
		int next_up = surf.Child(i + 1).PairI() - surf.Child(i).PairJ() - 1;
		sum_asymmetry += abs(prev_up - next_up);
		prev_up = next_up;
	}
	int next_up = surf.PairJ() - surf.Child(surf.NumChildren() - 1).PairJ() - 1;
	sum_asymmetry += abs(prev_up - next_up);
	sum_asymmetry += abs(next_up - first_up);
	return sum_asymmetry;
}

int librnary::SumAsymmetryUpperBound(unsigned num_nucleotides, unsigned min_hairpin_unpaired, unsigned branches) {
	// The worst configuration is always to clump them all together such that (_)(_)...(_)(_)
	// This works for circular and non-circular.
	return max(0, static_cast<int>(num_nucleotides - branches * (2 + min_hairpin_unpaired))) * 2;
}

string librnary::Stacking::ToString() const {
	stringstream ss;
	ss << "Stacking(" << i << ", " << j << ", ";
	if (t == MISMATCH)
		ss << "mismatch";
	else if (t == DANLGE5)
		ss << "5' dangle";
	else if (t == DANLGE3)
		ss << "3' dangle";
	else if (t == FLUSHCX)
		ss << "flush coax (" << k << "," << l << ")";
	else if (t == MMCX)
		ss << "mismatch-mediated coax (" << k << "," << l << ")";
	else
		ss << "none";
	ss << ")";
	return ss.str();
}

bool librnary::StackingFeature::IsBranch() const {
	return i != j;
}
bool librnary::StackingFeature::IsUnpaired() const {
	return !IsBranch();
}

bool librnary::IsBranch(const StackingFeature &f) {
	return f.IsBranch();
}

bool librnary::IsUnpaired(const StackingFeature &f) {
	return f.IsUnpaired();
}

vector<librnary::StackingFeature> librnary::MLStackingFeatures(const LoopRegion &loop) {
	int left_lim = loop.i;
	vector<StackingFeature> fts;
	for (const auto &e: loop.enclosed) {
		for (int i = left_lim + 1; i < e.i; ++i) {
			fts.push_back({i, i});
		}
		fts.push_back({e.i, e.j});
		left_lim = e.j;
	}
	for (int i = left_lim + 1; i < loop.j; ++i) {
		fts.push_back({i, i});
	}
	return fts;
}
void librnary::ExtractMultiLoops(std::vector<librnary::LoopRegion> &mls, const librnary::Surface &surf) {
	if (surf.NumChildren() >= 2 && !surf.IsExternalLoop()) {
		mls.push_back(surf);
	}
	for (const auto &cs : surf.Children()) {
		ExtractMultiLoops(mls, cs);
	}
}

int librnary::ExtractBranches(const librnary::LoopRegion &loop) {
	assert(loop.enclosed.size() >= 2);
	return (int) loop.enclosed.size() + 1;
}

std::vector<int> librnary::ExtractSegLengths(const librnary::LoopRegion &loop) {
	assert(loop.enclosed.size() >= 2);
	int last = loop.i;
	std::vector<int> segs;
	for (auto bp : loop.enclosed) {
		segs.push_back(bp.i - last - 1);
		last = bp.j;
	}
	segs.push_back(loop.j - last - 1);
	return segs;
}

int librnary::ExtractUnpaired(const librnary::LoopRegion &loop) {
	assert(loop.enclosed.size() >= 2);
	auto segs = ExtractSegLengths(loop);
	int count = 0;
	for (int s : segs)
		count += s;
	return count;
}

int librnary::ExtractSumAsymmetry(const librnary::LoopRegion &loop) {
	assert(loop.enclosed.size() >= 2);
	int sum_asymmetry = 0;
	int first_up = loop.enclosed[0].i - loop.i - 1;
	int prev_up = first_up;
	for (int i = 0; i < static_cast<int>(loop.enclosed.size()) - 1; ++i) {
		int next_up = loop.enclosed[i + 1].i - loop.enclosed[i].j - 1;
		sum_asymmetry += abs(prev_up - next_up);
		prev_up = next_up;
	}
	int next_up = loop.j - loop.enclosed.back().j - 1;
	sum_asymmetry += abs(prev_up - next_up);
	sum_asymmetry += abs(next_up - first_up);
	return sum_asymmetry;
}

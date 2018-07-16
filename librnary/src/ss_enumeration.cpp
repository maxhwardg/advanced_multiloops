//
// Created by max on 5/30/16.
//

#include <algorithm>

#include "ss_enumeration.hpp"

using namespace std;

vector<librnary::BondPair> librnary::BondPairs(const librnary::PrimeStructure &prim) {
	int sz = static_cast<int>(prim.size());
	vector<librnary::BondPair> bonds;
	for (int i = 0; i < sz; ++i) {
		for (int j = i + 1; j < sz; ++j) {
			if (librnary::ValidPair(prim[i], prim[j])) {
				bonds.emplace_back(i, j);
			}
		}
	}
	return bonds;
}

librnary::Matching librnary::RandomMatching(const librnary::PrimeStructure &prim, default_random_engine &re,
											unsigned max_trials, unsigned min_oneloop_unpaired) {
	auto bps = BondPairs(prim);
	shuffle(bps.begin(), bps.end(), re);
	auto match = EmptyMatching(static_cast<unsigned>(prim.size()));
	for (int i = 0; i < static_cast<int>(max_trials) && i < static_cast<int>(bps.size()); ++i) {
		auto bp = bps[i];
		if (match[bp.i] != bp.i || match[bp.j] != bp.j || (bp.j - bp.i - 1) < static_cast<int>(min_oneloop_unpaired)) {
			continue;
		}
		bool skip = false;
		for (int k = bp.i + 1; k < bp.j; ++k) {
			if (match[k] != k && (match[k] < bp.i || match[k] > bp.j)) {
				skip = true;
				break;
			}
		}
		if (!skip) {
			match[bp.i] = bp.j;
			match[bp.j] = bp.i;
		}
	}
	return match;
}

vector<vector<int>> librnary::ForwardBondTable(const PrimeStructure &bases) {
	vector<vector<int>> bonds(bases.size());
	for (unsigned i = 0; i < bases.size(); ++i)
		for (unsigned j = i + 1; j < bases.size(); ++j)
			if (ValidPair(bases[i], bases[j]))
				bonds[i].push_back((int) j);
	return bonds;
}

vector<vector<int>> librnary::BackwardBondTable(const PrimeStructure &bases) {
	vector<vector<int>> bonds(bases.size());
	for (unsigned i = 0; i < bases.size(); ++i)
		for (unsigned j = i + 1; j < bases.size(); ++j)
			if (ValidPair(bases[i], bases[j]))
				bonds[j].push_back((int) i);
	return bonds;
}

bool librnary::StructureEnumerator::MeetsConstraints(const PrimeStructure &bases, int i, int j) {
	return librnary::ValidPair(bases[i], bases[j]) && j - i - 1 >= min_hairpin;
}

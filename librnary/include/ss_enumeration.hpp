//
// Created by max on 5/30/16.
//

/**
 * This file contains functions and classes for enumerating secondary structure elements.
 */

#ifndef RNARK_SS_ENUMERATION_HPP
#define RNARK_SS_ENUMERATION_HPP

#include <stack>

#include "secondary_structure.hpp"
#include "models/nn_model.hpp"

namespace librnary {

/// This class does brute force enumeration of all secondary structures.
class StructureEnumerator {
	/// Minimum size of a hairpin loop in number of unpaired bases.
	int min_hairpin;
	/// The length of the RNA primary sequence being enumerated.
	int N{};
	/// Stack used during enumeration.
	std::stack<int> s;
	Matching match;

	/// Returns true if the bond i,j is allowable in a secondary structure.
	/// Allowable means that all the constraints for this enumerator are met.
	/// Also, it must be a valid RNA pairing.
	bool MeetsConstraints(const PrimeStructure &bases, int i, int j);

	template<typename Functor>
	void EnumerateHelper(const PrimeStructure &bases, Functor &f) {
		// Declaring these as local variables is a micro optimization, but it helps.
		int unmatched = static_cast<int>(s.size()), matched = static_cast<int>(match.size());
		// If there are more things on the stack than we can possibly match.
		// Then this is an impossible matching.
		if (unmatched > N - matched)
			return;

		if (matched == N) {
			f(match);
		} else {
			s.push((int) match.size());
			match.push_back(matched);
			EnumerateHelper(bases, f);
			match.pop_back();
			s.pop();
			if (unmatched && MeetsConstraints(bases, s.top(), matched)) {
				match[s.top()] = matched;
				match.push_back(s.top());
				s.pop();
				EnumerateHelper(bases, f);
				s.push(match[match.size() - 1]);
				match.pop_back();
			}
			match.push_back(matched);
			EnumerateHelper(bases, f);
			match.pop_back();
		}
	}

public:
	explicit StructureEnumerator(unsigned _min_hairpin)
		: min_hairpin(_min_hairpin) {}

	/// Apply some function f to a all secondary structures for a primary sequence "bases".
	template<typename Functor>
	void Enumerate(const PrimeStructure &bases, Functor &f) {
		s = std::stack<int>();
		match = Matching();
		N = (int) bases.size();
		EnumerateHelper(bases, f);
	}
};
std::vector<BondPair> BondPairs(const PrimeStructure &prim);
/**
 * Generates a random matching from a primary sequence. Attempts max_trials random valid bonds from the sequence.
 * Never generates pseudoknots.
 * Takes O(N^2+max_trials*N) time as it generates all valid bonds first. This ensures a good spread of random bonds.
 * Forces one-loops to have at least min_oneloop_unpaired unpaired nucleotides.
 */
Matching RandomMatching(const PrimeStructure &prim, std::default_random_engine &re, unsigned max_trials,
						unsigned min_oneloop_unpaired = NNModel::MIN_HAIRPIN_UNPAIRED);
std::vector<std::vector<int>> ForwardBondTable(const PrimeStructure &bases);
std::vector<std::vector<int>> BackwardBondTable(const PrimeStructure &bases);
}
#endif //RNARK_SS_ENUMERATION_HPP

//
// Created by max on 5/4/16.
//

/**
 * Contains types and functions for secondary structure.
 * We define secondary structure as the result of base pairings,
 * but excluding pseudoknots for such pairings.
 */

#ifndef RNARK_SECONDARY_STRUCTURE_HPP
#define RNARK_SECONDARY_STRUCTURE_HPP

#include <vector>
#include <memory>
#include <string>
#include <cassert>

#include "primary_structure.hpp"

namespace librnary {

/**
 * Unique IDs for RNA secondary structure loop types.
 */
enum LoopType {
	MULTI = 0,
	BULGE = 1,
	INTERNAL = 2,
	HAIRPIN = 3,
	EXTERNAL = 4,
	STACK = 5,
	NUM_LOOP_TYPES = 6,
};

/**
 * Matching represents a secondary structure.
 * The index i in the matching corresponds to the base i is bonded with.
 * In the case that i is unbonded, it will be paired with itself.
 * Another way of thinking about this that that a Matching is a permutation of integers [0, matching.size()).
 * This gives Matching some nice mathematical properties. It is closed under the "find paired" operator.
 * In addition, this sets up an involution: matching[matching[i]] == i.
 */
typedef std::vector<int> Matching;

/**
 * Creates and returns Matching of a particular size in which every base is unpaired.
 */
Matching EmptyMatching(unsigned size);

/// Represents the paired nucleotides i and j. Assumed to be 5' to 3' so i<j.
struct BondPair {
	// i and j are the nucleotides. Assumes i < j.
	int i{0}, j{0};
	BondPair(int _i, int _j)
		: i(_i), j(_j) {
		assert(i < j);
	}
	/// A BondPair constructed this way represents "no bond"
	BondPair() = default;
	/// Returns true if i,j < other.j.
	bool operator<(const BondPair &other) const;
};

/// Represents a stem/helix in an RNA structure.
struct Stem {
	/// The pair that closes all the other pairs in the stem. Has the lowest i and greatest j value.
	BondPair outer_pair{};
	/// Number of pairs in the stem.
	int num_pairs{0};
	Stem(BondPair bp, int np) : outer_pair(bp), num_pairs(np) {}
};

Matching DotBracketToMatching(const std::string &db);
std::string MatchingToDotBracket(const Matching &match);
/// Whether two bases are a valid canonical bond pairing.
bool ValidPair(Base a, Base b);
Matching BondPairsToMatching(const std::vector<BondPair> &bonds, unsigned sz);
std::vector<BondPair> MatchingToBondPairs(const Matching &match);
bool WatsonCrick(Base a, Base b);
bool ContainsPseudoknot(const Matching &match);

/// Returns true iff the pair (i,j) has to be a lonley pair by nature of its neighbours.
bool MustBeLonelyPair(const PrimeStructure &rna, int i, int j, int min_hairpin_unpaired);

/**
 * Finds all the stems in an RNA secondary structure. See the Stem struct.
 * @param match A matching representing a secondary structure.
 * @return A list of all the stems.
 */
std::vector<Stem> ExtractStems(const Matching &match);
}

#endif //RNARK_SECONDARY_STRUCTURE_HPP

//
// Created by max on 5/4/16.
//

#include <stack>
#include <algorithm>

#include "secondary_structure.hpp"

using namespace std;

librnary::Matching librnary::EmptyMatching(unsigned size) {
	Matching match(size);
	for (unsigned i = 0; i < size; ++i) {
		match[i] = i;
	}
	return match;
}

bool librnary::BondPair::operator<(const BondPair &other) const {
	return other.i > j;
}

librnary::Matching librnary::DotBracketToMatching(const string &db) {
	stack<int> s;
	auto match = EmptyMatching(static_cast<unsigned>(db.size()));

	for (int i = 0; i < static_cast<int>(db.size()); ++i) {
		if (db[i] == ')') {
			assert(!s.empty());
			match[i] = s.top();
			match[s.top()] = i;
			s.pop();
		} else if (db[i] == '(')
			s.push(i);
	}
	assert(s.empty());
	return match;
}

string librnary::MatchingToDotBracket(const librnary::Matching &match) {
	string db(match.size(), '.');
	for (int i = 0; i < static_cast<int>(match.size()); ++i) {
		if (match[i] > i)
			db[i] = '(';
		else if (match[i] < i)
			db[i] = ')';
	}
	return db;
}

bool librnary::ValidPair(librnary::Base a, librnary::Base b) {
	switch (a) {
		case A:
			return b == U;
		case U:
			switch (b) {
				case A:
				case G:
					return true;
				default:
					return false;
			}
		case G:
			switch (b) {
				case C:
				case U:
					return true;
				default:
					return false;
			}
		case C:
			return b == G;
		default:
			return false;
	}
}

bool librnary::WatsonCrick(librnary::Base a, librnary::Base b) {
	switch (a) {
		case A:
			return b == U;
		case U:
			return b == A;
		case G:
			return b == C;
		case C:
			return b == G;
		default:
			return false;
	}
}

librnary::Matching librnary::BondPairsToMatching(const vector<librnary::BondPair> &bonds, unsigned sz) {
	auto match = EmptyMatching(sz);
	for (auto &b : bonds) {
		match[b.i] = b.j;
		match[b.j] = b.i;
	}
	return match;
}

vector<librnary::BondPair> librnary::MatchingToBondPairs(const librnary::Matching &match) {
	vector<librnary::BondPair> bonds;
	for (int i = 0; i < static_cast<int>(match.size()); ++i)
		if (match[i] > i)
			bonds.push_back(BondPair(i, match[i]));
	return bonds;
}

bool librnary::ContainsPseudoknot(const librnary::Matching &match) {
	stack<int> open;
	for (int i = 0; i < static_cast<int>(match.size()); ++i) {
		if (match[i] > i) {
			open.push(i);
		} else if (match[i] < i) {
			assert(!open.empty());
			if (match[i] == open.top()) {
				open.pop();
			} else {
				return true;
			}
		}
	}
	assert(open.empty());
	return false;
}

bool librnary::MustBeLonelyPair(const PrimeStructure &rna, int i, int j, int min_hairpin_unpaired) {
	assert(i >= 0 && i < static_cast<int>(rna.size()) && j >= 0 && j < static_cast<int>(rna.size()));
	return !(i - 1 >= 0 && j + 1 < static_cast<int>(rna.size()) && ValidPair(rna[i - 1], rna[j + 1]))
		&& !(j - i - 3 >= min_hairpin_unpaired && ValidPair(rna[i + 1], rna[j - 1]));
}
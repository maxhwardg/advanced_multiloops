//
// Created by max on 5/4/16.
//

#include <stack>

#include "secondary_structure.hpp"

using namespace std;

namespace librnary {

Matching EmptyMatching(unsigned size) {
	Matching match(size);
	for (unsigned i = 0; i < size; ++i) {
		match[i] = i;
	}
	return match;
}

bool BondPair::operator<(const BondPair &other) const {
	return other.i > j;
}

Matching DotBracketToMatching(const string &db) {
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

string MatchingToDotBracket(const Matching &match) {
	string db(match.size(), '.');
	for (int i = 0; i < static_cast<int>(match.size()); ++i) {
		if (match[i] > i)
			db[i] = '(';
		else if (match[i] < i)
			db[i] = ')';
	}
	return db;
}

bool ValidPair(Base a, Base b) {
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

bool WatsonCrick(Base a, Base b) {
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

Matching BondPairsToMatching(const vector<BondPair> &bonds, unsigned sz) {
	auto match = EmptyMatching(sz);
	for (auto &b : bonds) {
		match[b.i] = b.j;
		match[b.j] = b.i;
	}
	return match;
}

vector<BondPair> MatchingToBondPairs(const Matching &match) {
	vector<BondPair> bonds;
	for (int i = 0; i < static_cast<int>(match.size()); ++i)
		if (match[i] > i)
			bonds.emplace_back(i, match[i]);
	return bonds;
}

bool ContainsPseudoknot(const Matching &match) {
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

bool MustBeLonelyPair(const PrimeStructure &rna, int i, int j, int min_hairpin_unpaired) {
	assert(i >= 0 && i < static_cast<int>(rna.size()) && j >= 0 && j < static_cast<int>(rna.size()));
	return !(i - 1 >= 0 && j + 1 < static_cast<int>(rna.size()) && ValidPair(rna[i - 1], rna[j + 1]))
		&& !(j - i - 3 >= min_hairpin_unpaired && ValidPair(rna[i + 1], rna[j - 1]));
}

std::vector<Stem> ExtractStems(const Matching &match) {
	std::vector<Stem> stems;
	std::vector<int> marked(match.size());
	for (int i = 0; i < static_cast<int>(match.size()); ++i) {
		if (marked[i])
			continue;
		if (match[i] == i)
			continue;
		marked[i] = marked[match[i]] = 1;
		int ip = i, jp = match[i], sz = 1;
		while (ip + 1 < jp - 1 && match[ip + 1] == jp - 1) {
			++ip;
			--jp;
			++sz;
			marked[ip] = marked[jp] = 1;
		}
		stems.emplace_back(BondPair(i, match[i]), sz);
	}
	return stems;
}

}

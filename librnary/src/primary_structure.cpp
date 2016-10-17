//
// Created by max on 5/4/16.
//

#include "primary_structure.hpp"

using namespace std;

librnary::Base librnary::CharToBase(char c) {
	// Note that we don't toupper(c) as a lower case
	// indicates a non-pairing nucleotide.
	switch (c) {
		case 'A':
			return A;
		case 'T':
			return U;
		case 'U':
			return U;
		case 'G':
			return G;
		case 'C':
			return C;
		default:
			break;
	}
	return NUMBASES; // Invalid base.
}

char librnary::BaseToChar(librnary::Base b) {
	switch (b) {
		case A:
			return 'A';
		case U:
			return 'U';
		case G:
			return 'G';
		case C:
			return 'C';
		case NUMBASES:
			return 'X'; // Special case for invalid base.
	}
	return 'X'; // this is not possible
}

librnary::PrimeStructure librnary::StringToPrimary(const string &seq) {
	librnary::PrimeStructure bases;
	bases.reserve(seq.size());
	for (char c : seq)
		bases.push_back(CharToBase(c));
	return bases;
}

std::string librnary::PrimaryToString(const PrimeStructure &primary) {
	string str;
	for (auto b : primary)
		str.push_back(BaseToChar(b));
	return str;
}
librnary::PrimeStructure librnary::RandomPrimary(default_random_engine &re, unsigned length) {
	PrimeStructure ps(length);
	vector<Base> bases = {A, U, G, C};
	uniform_int_distribution<int> dist(0, static_cast<int>(bases.size() - 1));
	for (Base &b : ps) {
		b = bases[dist(re)];
	}
	return ps;
}



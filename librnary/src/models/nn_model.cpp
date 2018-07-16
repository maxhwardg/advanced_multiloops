//
// Created by max on 5/31/16.
//

#include "models/nn_model.hpp"

using namespace std;

librnary::energy_t librnary::NNModel::OneLoop(int i, int j) const {
	assert(i < j);
	assert(ValidPair(rna[i], rna[j]));
	if (j - i - 1 < MIN_HAIRPIN_UNPAIRED)
		return MaxMFE();
	return erg3(i + 1, j + 1, struc.get(), dt.get(), 0);
}

librnary::energy_t librnary::NNModel::TwoLoop(int i, int k, int l, int j) const {
	assert(i < j && k < l && i < k && l < j);
	if (k == i + 1 && l == j - 1)
		return erg1(i + 1, j + 1, k + 1, l + 1, struc.get(), dt.get());
	return erg2(i + 1, j + 1, k + 1, l + 1, struc.get(), dt.get(), 0, 0);
}

librnary::energy_t librnary::NNModel::Branch(int i, int j) const {
	return penalty(i + 1, j + 1, struc.get(), dt.get());
}


librnary::energy_t librnary::NNModel::FlushCoax(int i, int j, int k, int l) const {
	assert(i < j && k < l && (k == j + 1 || k == i + 1 || l == j - 1));
	if (i < k && l < j) // i,j close a multi-loop.
		return erg1(i + 1, j + 1, k + 1, l + 1, struc.get(), dt.get());
	return erg1(j + 1, i + 1, k + 1, l + 1, struc.get(), dt.get());
}

librnary::energy_t librnary::NNModel::MismatchCoax(int i, int j, int k, int l) const {
	assert(i < j && k < l && (abs(i - k) == 2 || abs(i - l) == 2 || abs(j - k) == 2 || abs(j - l) == 2));
	// i,j is mismatched stacked with k,l
	// Mismatch is off i,j
	i += 1;
	j += 1;
	k += 1;
	l += 1;
	if (i < k && l < j) { // i,j closes a multiloop
		if (k == i + 2) // (.(_)_.)
			return ergcoaxinterbases1(j, i, k, l, struc.get(), dt.get());
		// (._(_).)
		return ergcoaxinterbases2(k, l, j, i, struc.get(), dt.get());
	} else if (k < i && j < l) { // k,l closes a multiloop
		if (i == k + 2) // (.(_)._)
			return ergcoaxinterbases2(l, k, i, j, struc.get(), dt.get());
		// (_.(_).)
		return ergcoaxinterbases1(i, j, l, k, struc.get(), dt.get());
	} else if (j < k) { // .(_).(_)
		return ergcoaxinterbases1(i, j, k, l, struc.get(), dt.get());
	} else { // (_).(_).
		return ergcoaxinterbases2(k, l, i, j, struc.get(), dt.get());
	}

}

librnary::energy_t librnary::NNModel::FiveDangle(int i, int j) const {
	assert(i < j && i > 0);
	return erg4(j + 1, i + 1, i, 2, struc.get(), dt.get(), false);
}

librnary::energy_t librnary::NNModel::ClosingFiveDangle(int i, int j) const {
	assert(i < j);
	return erg4(i + 1, j + 1, j, 2, struc.get(), dt.get(), false);
}

librnary::energy_t librnary::NNModel::ThreeDangle(int i, int j) const {
	assert(i < j && j + 1 < static_cast<int>(rna.size()));
	return erg4(j + 1, i + 1, j + 2, 1, struc.get(), dt.get(), false);
}

librnary::energy_t librnary::NNModel::ClosingThreeDangle(int i, int j) const {
	assert(i < j);
	return erg4(i + 1, j + 1, i + 2, 1, struc.get(), dt.get(), false);
}

librnary::energy_t librnary::NNModel::Mismatch(int i, int j) const {
	assert(i < j);
	return dt->tstkm[struc->numseq[j + 1]][struc->numseq[i + 1]][struc->numseq[j + 2]][struc->numseq[i]];
}

librnary::energy_t librnary::NNModel::ClosingMismatch(int i, int j) const {
	assert(i < j);
	return dt->tstkm[struc->numseq[i + 1]][struc->numseq[j + 1]][struc->numseq[i + 2]][struc->numseq[j]];
}

void librnary::NNModel::SetRNA(const librnary::PrimeStructure &primary) {
	this->rna = primary;
	this->struc = librnary::LoadStructure(rna);
}

librnary::PrimeStructure librnary::NNModel::RNA() const {
	return this->rna;
}

librnary::energy_t librnary::NNModel::MaxMFE() const {
	return numeric_limits<librnary::energy_t>::max() / 3;
}
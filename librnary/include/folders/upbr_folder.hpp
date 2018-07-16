//
// Created by max on 6/22/16.
//

#ifndef RNARK_UPBR_FOLDER_HPP
#define RNARK_UPBR_FOLDER_HPP

#include "vector_types.hpp"
#include "models/nn_unpaired_model.hpp"

#include <stack>

namespace librnary {
/**
 * A folding algorithm that can incorporate multi-loop functions depending on #unpaired and #branches.
 */
class UpBrFolder {

	// P stands for "Paired" and is the same as the V table.
	// P[i][j]: MFE of structures closed by the bond i,j.
	VVE P;
	// The "Multiloop" table.
	// ML[b][up][i][j]: MFE of multiloop fragments from i to j having >= b branches and up unpaired nucleotides.
	_4DVE ML;
	// The external loop table.
	// E[i]: The MFE external loop from 0 to i.
	VE E;

	VVE CxFl;
	VVE CxMM;

	PrimeStructure rna;

	StrainedUnpairedModel em;

	enum Table {
		ET, MLT, PT, CxFlT, CxMMT
	};

	// Trace state during traceback. Represents a location in a table.
	struct TState {
		Table t;
		int b, up, i, j;

		TState(int _i)
			: t(ET), i(_i) {};

		TState(int _i, int _j)
			: t(PT), i(_i), j(_j) {};

		TState(Table _t, int _i, int _j)
			: t(_t), i(_i), j(_j) {};

		TState(int _b, int _up, int _i, int _j)
			: t(MLT), b(_b), up(_up), i(_i), j(_j) {};
	};

	// Traceback through a place in the E table. Assumes this place/state is at the top of s.
	// Updates s accordingly with the decomposed states.
	void TraceE(std::stack<TState> &s);

	// Traceback through a place in the B table. Assumes this place/state is at the top of s.
	// Updates s accordingly with the decomposed states.
	void TraceP(std::stack<TState> &s);

	void TraceCxFl(std::stack<TState> &s);

	void TraceCxMM(std::stack<TState> &s);

	// Traceback through a place in the ML table. Assumes this place/state is at the top of s.
	// Updates s accordingly with the decomposed states.
	void TraceML(std::stack<TState> &s);

	/*
	 * The optimal sub-surface score of the structure closed by (i,j). Designed for external-loop branches.
	 * Accounts for AU/GU penalty.
	 */
	energy_t SSScore(int i, int j);

	// The maximum number of unpaired nucleotides in a multiloop.
	int max_multi_unpaired = std::numeric_limits<int>::max() / 3;
	int max_multi_branches = std::numeric_limits<int>::max() / 3;
	// The maximum number of unpaired in a bulge or internal loop.
	int max_twoloop_unpaired = std::numeric_limits<int>::max() / 3;

	// This flag toggles whether lonely pairs are allowed.
	bool lonely_pairs = true;

public:

	void SetMaxBranches(unsigned v);

	int MaxBranches() const;

	void SetLonelyPairs(bool v);

	bool LonelyPairs() const;

	/// Set the max unpaired nucleotides in a multiloop.
	void SetMaxMulti(unsigned max_up);

	/// Get the max unpaired nucleotides in a multiloop.
	int MaxMulti() const;

	/// Set the max unpaired nucleotides in a bulge/internal loop.
	void SetMaxTwoLoop(unsigned max_up);

	/// Get the max unpaired nucleotides in a bulge/internal loop.
	int MaxTwoLoop() const;

	UpBrFolder(const StrainedUnpairedModel &_em)
		: em(_em) {}

	/// Fill the DP tables and return the MFE value.
	virtual energy_t Fold(const PrimeStructure &rna);

	/// Trace back through the DP tables and produce a dot bracket representation.
	virtual Matching Traceback();

	VVE GetP() const;

	VE GetE() const;

	// Get a copy of the energy model.
	StrainedUnpairedModel GetEM() const;

};
}

#endif //RNARK_UPBR_FOLDER_HPP

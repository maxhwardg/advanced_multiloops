//
// Created by max on 6/19/16.
//

#ifndef RNARK_NN_UNPAIRED_FOLDER_HPP
#define RNARK_NN_UNPAIRED_FOLDER_HPP

#include "models/nn_unpaired_model.hpp"
#include "vector_types.hpp"

#include <stack>

namespace librnary {
/**
 * This class implements a dynamic programming algorithm for computing a MFE folding on an RNA sequence under the
 * nearest neighbour model with an arbitrary function on unpaired nucleotides in a multi-loop.
 * The complexity of folding is O(N^4).
 */
class NNUnpairedFolder {
	/// P stands for "Paired" and is the same as the V table.
	/// P[i][j]: MFE of structures closed by the bond i,j.
	VVE P;

	/// Rather bafflingly, nested vectors seem faster than flat arrays here.
	/// Declaring a large flat array seems very exnpensive. Plus the access pattern seems to mess with the pre-fetcher.
	/// The "Multiloop" table.
	/// ML[b][up][i][j]: MFE of multiloop fragments from i to j having >= b branches and up unpaired nucleotides.
	_4DVE ML;
	/// The external loop table.
	/// E[i]: The MFE external loop from 0 to i.
	VE E;

	VVE CxFl;
	VVE CxMM;

	PrimeStructure rna;

	NNUnpairedModel em;

	enum Table {
		ET, MLT, PT, CxFlT, CxMMT
	};

	/// Trace state during traceback. Represents a location in a table.
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

	/// Traceback through a place in the E table. Assumes this place/state is at the top of s.
	/// Updates s accordingly with the decomposed states.
	void TraceE(std::stack<TState> &s);

	/// Traceback through a place in the P table. Assumes this place/state is at the top of s.
	/// Updates s accordingly with the decomposed states.
	void TraceP(std::stack<TState> &s);

	void TraceCxFl(std::stack<TState> &s);

	void TraceCxMM(std::stack<TState> &s);

	/// Traceback through a place in the ML table. Assumes this place/state is at the top of s.
	/// Updates s accordingly with the decomposed states.
	void TraceML(std::stack<TState> &s);

	/**
	 * The optimal sub-surface score of the structure closed by (i,j). Designed for external-loop branches.
	 * Accounts for AU/GU penalty.
	 */
	energy_t SSScore(int i, int j) const;
	/**
	 * The optimal sub-surface score of the structure closed by (i,j) in a multi-loop.
	 * Accounts for AU/GU penalty, and multi-loop branch cost.
	 */
	energy_t MLSSScore(int i, int j) const;

	/// The maximum number of unpaired nucleotides in a multiloop.
	int max_multi_unpaired = std::numeric_limits<int>::max() / 3;
	/// The maximum number of unpaired in a bulge or internal loop.
	int max_twoloop_unpaired = std::numeric_limits<int>::max() / 3;

	/// This flag toggles whether lonely pairs are allowed.
	bool lonely_pairs = true;

	/// This flag toggles whether stacking interactions (dangles, terminal mismatch, and coaxial stacking) are used.
	bool stacking = true;

public:

	void SetStacking(bool v);

	bool Stacking() const;

	void SetLonelyPairs(bool v);

	bool LonelyPairs() const;

	/// Set the max unpaired nucleotides in a multiloop.
	void SetMaxMulti(unsigned max_up);

	/// Get the max unpaired nucleotides in a multiloop.
	int MaxMultiUnpaired() const;

	/// Set the max unpaired nucleotides in a bulge/internal loop.
	void SetMaxTwoLoop(unsigned max_up);

	/// Get the max unpaired nucleotides in a bulge/internal loop.
	int MaxTwoLoop() const;

	/**
	 * Note that this method may make the previous call to fold invalid, as it will still assume the old model.
	 * @param _em Sets the model to use internally to this.
	 */
	void SetModel(const NNUnpairedModel &_em);

	NNUnpairedFolder(const NNUnpairedModel &_em)
		: em(_em) {}

	/// Fill the DP tables and return the MFE value.
	energy_t Fold(const PrimeStructure &rna);

	/// Trace back through the DP tables and produce a MFE secondary structure.
	Matching Traceback();

	VVE GetP() const;

	VE GetE() const;

	/// Get a copy of the energy model.
	NNUnpairedModel GetEM() const;

};
}

#endif //RNARK_NN_UNPAIRED_FOLDER_HPP

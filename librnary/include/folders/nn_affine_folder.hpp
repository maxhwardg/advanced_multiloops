//
// Created by max on 7/24/17.
//

#ifndef RNARK_NN_AFFINE_FOLDER_HPP
#define RNARK_NN_AFFINE_FOLDER_HPP

#include <models/nn_affine_model.hpp>
#include <energy.hpp>
#include <vector_types.hpp>
#include <stack>

namespace librnary {
/**
 * This class does MFE RNA folding using the typical Turner-Mathews nearest neighbour rules. However, it also
 * incorporates a particular 2017 Honying model of single nucleotide bulges. This model uses the most and least
 * stable stem for a single nucleotide bulge.
 */
class NNAffineFolder {
protected:
	PrimeStructure rna;
	NNAffineModel em;

	/**
	 * The 'Paired' table. P[i][j] is the optimal substructure closed by a pair i,j.
	 */
	VVE P;

	/**
	 * The 'Multi-Loop' table. ML[b][i][j] is the optimal part of the multi-loop that definitely has
	 * at least b branches.
	 */
	VVVE ML;

	/**
	 * The 'Coaxial stack' table. Cx[i][j] is the optimal multi-loop coaxial stack such that one branch starts at i,
	 * and the other branch ends at j. Note that it is assumed that these branches are in a multi-loop. This table will
	 * count the unpaired cost.
	 */
	VVE Cx;

	/**
	 * The 'External loop' table. E[i] is the optimal external loop fragment 0..i.
	 */
	VE E;

	/**
	 * Represents a particular table during a traceback.
	 */
	enum Table {
		ET, PT, MLT, CxT
	};

	/**
	 * Trace state during traceback. Represents a location in a table.
	 */
	struct TState {
		Table t;
		int extra, i, j;

		TState(int _i)
			: t(ET), i(_i) {};

		TState(Table _t, int _i, int _j)
			: t(_t), i(_i), j(_j) {};

		TState(Table _t, int _extra, int _i, int _j)
			: t(_t), extra(_extra), i(_i), j(_j) {};
	};

	/**
	 * Traceback through a place in the E table. Assumes this place/state is at the top of s.
	 * Updates s accordingly with the decomposed states.
	 * @param s Current open states in the traceback in stack form.
	 */
	void TraceE(std::stack<TState> &s);

	/**
	 * Traceback through a place in the P table. Assumes this place/state is at the top of s.
	 * Updates s accordingly with the decomposed states.
	 * @param s Current open states in the traceback in stack form.
	 */
	void TraceP(std::stack<TState> &s);

	/**
	 * Traceback through a place in the Cx table. Assumes this place/state is at the top of s.
	 * Updates s accordingly with the decomposed states.
	 * @param s Current open states in the traceback in stack form.
	 */
	void TraceCx(std::stack<TState> &s);

	/**
	 * Traceback through a place in the ML table. Assumes this place/state is at the top of s.
	 * Updates s accordingly with the decomposed states.
	 * @param s Current open states in the traceback in stack form.
	 */
	void TraceML(std::stack<TState> &s);

	/// This flag toggles whether lonely pairs are allowed.
	bool lonely_pairs = true;


	/// This flag toggles whether stacking interactions (dangles, terminal mismatch, and coaxial stacking) are used.
	bool stacking = true;

	/// The maximum number of unpaired in a bulge or internal loop.
	int max_twoloop_unpaired = std::numeric_limits<int>::max() / 3;

	/**
	 * The optimal sub-surface score of the structure closed by (i,j). Designed for external-loop branches.
	 * Accounts for AU/GU penalty.
	 */
	energy_t SSScore(int i, int j) const;

	/**
	 * The optimal sub-surface score of the structure closed by (i,j) in a multi-loop.
	 * Accounts for AU/GU penalty, and multi-loop branch cost.
	 */
	virtual energy_t MLSSScore(int i, int j) const;

	/**
	 * Computes Multi-loop closure free energy change. Also includes the cost of the closing branch.
	 * @param i 5' nucleotide of closing pair for multi-loop.
	 * @param j 3' nucleotide of closing pair for multi-loop.
	 * @return The free energy change of multi-loop closure.
	 */
	virtual energy_t MLClosingBranchScore(int i, int j) const;

public:

	/// Set the max unpaired nucleotides in a bulge/internal loop.
	void SetMaxTwoLoop(unsigned max_up);

	/// Get the max unpaired nucleotides in a bulge/internal loop.
	int MaxTwoLoop() const;

	void SetStacking(bool v);

	bool Stacking() const;

	void SetLonelyPairs(bool v);

	bool LonelyPairs() const;

	energy_t Fold(const PrimeStructure &_rna);

	/**
	 * Note that this method may make the previous call to fold invalid, as it will still assume the old model.
	 * @param _em Sets the model to use internally to this.
	 */
	void SetModel(const NNAffineModel &_em);

	/**
	 * Trace back through the DP tables and produce a MFE secondary structure.
	 * @return The traced secondary structure.
	 */
	Matching Traceback();

	NNAffineFolder(const NNAffineModel &_em)
		: em(_em) {}

};

}

#endif //RNARK_NN_AFFINE_FOLDER_HPP
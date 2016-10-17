//
// Created by max on 7/9/16.
//

#ifndef RNARK_ASYMMETRY_FOLDER_HPP
#define RNARK_ASYMMETRY_FOLDER_HPP

#include "asymmetry_model.hpp"
#include "vector_types.hpp"
#include "multi_array.hpp"

#include <stack>

namespace librnary {
class AsymmetryFolder {
	AsymmetryModel em;

	PrimeStructure rna;

	// DP Tables.

	/// The External loop table.
	VE E;
	/// The paired table.
	VVE P;
	/// The Coaxial Flush table.
	VVE CxFl;
	/// The Coaxial Mismatch 5' unpaired table.
	VVE CxMM5;
	/// The Coaxial Mismatch 3' unpaired table.
	VVE CxMM3;

	/// These arrays are hybrid flat/nested vector to optimize speed. This is empirically
	/// faster than all flat, or all nested vectors, or other hybrids.
	/// The Multi-Loop Unpaired table
	Array2D<_4DVE> ML_Up;
	/// The Multi-Loop Branch table.
	Array2D<_4DVE> ML_Br;


	enum Table {
		ET, PT, ML_UpT, ML_BrT, CxFlT, CxMM5T, CxMM3T
	};

	/// Trace state during traceback. Represents a location in a table.
	struct TState {
		Table t;
		int br, used_mask, up, upr, i, j;

		TState(int _i)
			: t(ET), i(_i) {};

		TState(int _i, int _j)
			: t(PT), i(_i), j(_j) {};

		TState(Table _t, int _i, int _j)
			: t(_t), i(_i), j(_j) {
			assert(t == CxFlT || t == CxMM5T || t == CxMM3T || t == PT);
		};

		TState(Table _t, int _br, int _used_mask, int _up, int _upr, int _i, int _j)
			: t(_t), br(_br), used_mask(_used_mask), up(_up), upr(_upr), i(_i), j(_j) {
			assert(t == ML_BrT || t == ML_UpT);
		};
	};

	int max_unpaired_gap = std::numeric_limits<int>::max() / 3;
	int max_twoloop_unpaired = std::numeric_limits<int>::max() / 3;
	bool lonely_pairs = true;
	bool stacking = true;

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

	energy_t AsymScore(int a, int b) const;

	/**
	 * @return Default minimum number of branches to force the DP to place in the internal part of a multi-loop.
	 */
	int DefaultRequiredMultiInternalBranches() const;

	/**
	 * @return Longest allowed gap of consecutive unpaired nucleotides in a multi-loop.
	 */
	int UnpairedGapLimit() const;

	// These traceback helper functions assume their place/state is at the top of s.
	// They update s accordingly with the decomposed states.

	void TraceE(std::stack<TState> &s);
	void TraceP(std::stack<TState> &s);
	void TraceCxFl(std::stack<TState> &s);
	void TraceCxMM3(std::stack<TState> &s);
	void TraceCxMM5(std::stack<TState> &s);
	void TraceMLUp(std::stack<TState> &s);
	void TraceMLBr(std::stack<TState> &s);

	void Relax(Table parent, energy_t &best, std::vector<TState> &best_decomp,
			   const std::vector<TState> &decomp, energy_t aux_e) const;

public:
	energy_t Fold(const PrimeStructure &rna);
	/// Trace-back through the DP tables and produce the secondary structure.
	Matching Traceback();
	AsymmetryFolder(const AsymmetryModel &_em)
		: em(_em) {}

	/// Gets the maximum allowed gap of consecutive unpaired nucleotides allowed in a multi-loop.
	int UnpairedGap() const;

	/// Sets the maximum allowed gap of consecutive unpaired nucleotides allowed in a multi-loop.
	void SetUnpairedGap(unsigned v);

	void SetStacking(bool v);

	bool Stacking() const;

	void SetLonelyPairs(bool v);

	bool LonelyPairs() const;

	/// Set the max unpaired nucleotides in a bulge/internal loop.
	void SetMaxTwoLoop(unsigned max_up);

	/// Get the max unpaired nucleotides in a bulge/internal loop.
	int MaxTwoLoop() const;


	VVE GetP() const;
};
}

#endif //RNARK_ASYMMETRY_FOLDER_HPP

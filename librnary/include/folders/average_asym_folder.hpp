//
// Created by max on 8/10/16.
//

#ifndef RNARK_AVERAGE_ASYM_FOLDER_HPP_HPP
#define RNARK_AVERAGE_ASYM_FOLDER_HPP_HPP

#include "vector_types.hpp"
#include "models/average_asym_model.hpp"

#include <stack>

namespace librnary {
/**
 * An MFE folding algorithm that incorporates the average asymmetry model.
 * This mode comes from "Experimentally derived nearest-neighbor parameters
 * for the stability of RNA three-and four-way multibranch loops" by Mathews & Turner in 2002.
 */
class AverageAsymmetryFolder {
	AverageAsymmetryModel em;

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
	/// The Multi-Loop Unpaired table
	/// ML_Up[bs][br][sum_asym][up_left][up_right][i][j] =
	/// The optimal multi-loop fragment between i and j given that 0 or more consecutive unpaired start at i.
	/// Also, the fragment has exactly br branches, and has exactly sum_asym asymmetry.
	/// Also, has up_left unpaired nucleotides for the preceding branch.
	/// Also has up_right unpaired nucleotides between the closing branch and the rightmost internal branch.
	/// Finally, the bs bitset represents whether certain unpaired nucleotides were used.
	_7DVE ML_Up;
	/// The Multi-Loop Branch table.
	/// Similar to ML_Up but assumes a branch starts at i.
	_7DVE ML_Br;


	enum Table {
		ET, PT, ML_UpT, ML_BrT, CxFlT, CxMM5T, CxMM3T
	};

	/// Trace state during traceback. Represents a location in a table.
	struct TState {
		Table t;
		int used_mask, br, sum_asym, up, upr, i, j;

		TState(int _i)
			: t(ET), i(_i) {};

		TState(int _i, int _j)
			: t(PT), i(_i), j(_j) {};

		TState(Table _t, int _i, int _j)
			: t(_t), i(_i), j(_j) {
			assert(t == CxFlT || t == CxMM5T || t == CxMM3T || t == PT);
		};

		TState(Table _t, int _used_mask, int _br, int _sum_asym, int _up, int _upr, int _i, int _j)
			: t(_t), used_mask(_used_mask), br(_br), sum_asym(_sum_asym), up(_up), upr(_upr), i(_i), j(_j) {
			assert(t == ML_BrT || t == ML_UpT);
		};
	};

	/// Max consecutive unpaired nucleotides allowed in a multi-loop internal fragment.
	int max_unpaired_gap = std::numeric_limits<int>::max() / 3;
	int max_twoloop_unpaired = std::numeric_limits<int>::max() / 3;
	/// Max branches allowed in a multi-loop internal fragment.
	int max_ml_branches = std::numeric_limits<int>::max() / 3;
	/// Maximum sum asymmetry a multi-loop internal fragment can have.
	int max_nonclosing_ml_sum_asym_ = std::numeric_limits<int>::max() / 3;
	bool lonely_pairs = true;
	bool stacking = true;

	/**
	 * @return The upper bound on number of branches in a multi-loop.
	 */
	int BranchesUB() const;

	/**
	 * @return The upper bound on the sum of asymmetry for all but the closing branch in any multi-loop.
	 */
	int NonClosingSumAsymmetryUB() const;

	/**
	 * @return The upper bound on the length of a gap of consecutive unpaired nucleotides in a multi-loop.
	 */
	int UnpairedGapUB() const;


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
	AverageAsymmetryFolder(const AverageAsymmetryModel &_em)
		: em(_em) {}

	/// Gets the max number of consecutive unpaired nucleotides for the internal part of a multi-loop.
	int UnpairedGap() const;

	/// Sets the max number of consecutive unpaired nucleotides for the internal part of a multi-loop.
	void SetUnpairedGap(unsigned v);

	void SetStacking(bool v);

	bool Stacking() const;

	void SetLonelyPairs(bool v);

	bool LonelyPairs() const;

	/// Set the max unpaired nucleotides in a bulge/internal loop.
	void SetMaxTwoLoop(unsigned max_up);

	/// Get the max unpaired nucleotides in a bulge/internal loop.
	int MaxTwoLoop() const;

	/// Gets the max number of multi-loop branches for the internal part of a multi-loop.
	int MaxMLBranches() const;

	/// Sets the max number of multi-loop branches for the internal part of a multi-loop.
	void SetMaxMLBranches(int v);

	/**
	 * @return The lower bound in sum asymmetry for multi-loops.
	 */
	int MaxMLNonClosingAsym() const;
	/**
	 * @param v The new lower bound for sum asymmetry for multi-loops.
	 */
	void SetMaxMLNonClosingAsym(int v);


	VVE GetP() const;
};
}

#endif //RNARK_AVERAGE_ASYM_FOLDER_HPP_HPP

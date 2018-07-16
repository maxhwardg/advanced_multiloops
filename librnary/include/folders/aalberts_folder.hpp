//
// Created by max on 6/28/16.
//

#ifndef RNARK_AALBERTS_FOLDER_HPP
#define RNARK_AALBERTS_FOLDER_HPP

#include <stack>

#include "primary_structure.hpp"
#include "models/aalberts_model.hpp"
#include "vector_types.hpp"

namespace librnary {

/**
 * A folding algorithm that can use Aalberts-Nandagopal like multi-loop models.
 */
class AalbertsFolder {

	/// P stands for "Paired" and is the same as the V table.
	/// P[i][j]: MFE of structures closed by the bond i,j.
	VVE P;
	/// The "Multiloop" table.
	_5DVE ML;
	/// The external loop table.
	/// E[i]: The MFE external loop from 0 to i.
	VE E;
	/// The coaxial stacking table.
	VVE Cx;

	PrimeStructure rna;

	AalbertsModel em;

	enum Table {
		ET, MLT, PT, CxT
	};

	/// Trace state during traceback. Represents a location in a table.
	struct TState {
		Table t;
		int br, b, a, i, j;

		TState(int _i)
			: t(ET), i(_i) {};

		TState(int _i, int _j)
			: t(PT), i(_i), j(_j) {};

		TState(Table _t, int _i, int _j)
			: t(_t), i(_i), j(_j) {};

		TState(int _br, int _b, int _a, int _i, int _j)
			: t(MLT), br(_br), b(_b), a(_a), i(_i), j(_j) {};
	};

	/// Traceback through a place in the E table. Assumes this place/state is at the top of s.
	/// Updates s accordingly with the decomposed states.
	void TraceE(std::stack<TState> &s);

	/// Traceback through a place in the B table. Assumes this place/state is at the top of s.
	/// Updates s accordingly with the decomposed states.
	void TraceP(std::stack<TState> &s);

	void TraceCx(std::stack<TState> &s);

	/// Traceback through a place in the ML table. Assumes this place/state is at the top of s.
	/// Updates s accordingly with the decomposed states.
	void TraceML(std::stack<TState> &s);

	/**
	 * The optimal sub-surface score of the structure closed by (i,j). Designed for external-loop branches.
	 * Accounts for AU/GU penalty.
	 */
	energy_t SSScore(int i, int j);

	/// The maximum number of multi-loop a-length and b-length segments.
	int max_multi_alength = std::numeric_limits<int>::max() / 3;
	int max_multi_blength = std::numeric_limits<int>::max() / 3;
	/// The maximum number of unpaired in a bulge or internal loop.
	int max_twoloop_unpaired = std::numeric_limits<int>::max() / 3;

	/// Maximum number of B-Length segments in a range of n nucleotides.
	int MaxBLengthSegs(int n) const {
		return std::min(n / 2, max_multi_blength);
	}

	/// Maximum number of A-Length segments in a range of n nucleotides.
	int MaxALengthSegs(int n) const {
		return std::min(n, max_multi_alength);
	}

	/// This flag toggles whether lonely pairs are allowed.
	bool lonely_pairs = true;

public:

	/// Set the max length-b segments in the internal part of a multi-loop.
	void SetMaxBLength(unsigned v);

	/// Get the max length-b segments in the internal part of a multi-loop.
	int MaxBLength() const;

	void SetLonelyPairs(bool v);

	bool LonelyPairs() const;

	/// Set the max length-a segments in the internal part of a multi-loop.
	void SetMaxALength(unsigned max_up);

	/// Get the max length-a segments in the internal part of a multi-loop.
	int MaxALength() const;

	/// Set the max unpaired nucleotides in a bulge/internal loop.
	void SetMaxTwoLoop(unsigned max_up);

	/// Get the max unpaired nucleotides in a bulge/internal loop.
	int MaxTwoLoop() const;

	AalbertsFolder(const AalbertsModel &_em)
		: em(_em) {}

	/// Fill the DP tables and return the MFE value.
	virtual energy_t Fold(const PrimeStructure &rna);

	/// Trace back through the DP tables and produce a dot bracket representation.
	virtual Matching Traceback();

	VVE GetP() const;

	VE GetE() const;

	/// Get a copy of the energy model.
	AalbertsModel GetEM() const;

	void SetModel(const AalbertsModel &em);

};
}

#endif //RNARK_AALBERTS_FOLDER_HPP

//
// Created by max on 7/9/16.
//
#include "asymmetry_folder.hpp"

using namespace std;

librnary::energy_t librnary::AsymmetryFolder::SSScore(int i, int j) const {
	return em.Branch(i, j) + P[i][j];
}

librnary::energy_t librnary::AsymmetryFolder::MLSSScore(int i, int j) const {
	return SSScore(i, j) + em.MLBranchCost();
}

librnary::energy_t librnary::AsymmetryFolder::AsymScore(int a, int b) const {
	return em.MLAsymmetryCost() * abs(a - b);
}

int librnary::AsymmetryFolder::DefaultRequiredMultiInternalBranches() const {
	return 2;
}

int librnary::AsymmetryFolder::UnpairedGapLimit() const {
	// Since the largest gap in a multi-loop must look like this:
	// ((...)(...)_)
	// This admits a tight limit on the largest unpaired gap.
	return min(max_unpaired_gap, max(0, int(rna.size() - 2 - 2 * (2 + em.MIN_HAIRPIN_UNPAIRED))));
}


int librnary::AsymmetryFolder::UnpairedGap() const {
	return max_unpaired_gap;
}
void librnary::AsymmetryFolder::SetUnpairedGap(unsigned v) {
	max_unpaired_gap = v;
}


void librnary::AsymmetryFolder::SetMaxTwoLoop(unsigned max_up) {
	max_twoloop_unpaired = max_up;
}

int librnary::AsymmetryFolder::MaxTwoLoop() const {
	return max_twoloop_unpaired;
}

void librnary::AsymmetryFolder::SetLonelyPairs(bool v) {
	lonely_pairs = v;
}

bool librnary::AsymmetryFolder::LonelyPairs() const {
	return lonely_pairs;
}

void librnary::AsymmetryFolder::SetStacking(bool v) {
	stacking = v;
}

bool librnary::AsymmetryFolder::Stacking() const {
	return stacking;
}


librnary::VVE librnary::AsymmetryFolder::GetP() const {
	return P;
}

void librnary::AsymmetryFolder::Relax(Table parent, int &best, vector<TState> &best_decomp,
									  const vector<TState> &decomp, int aux_e) const {
	int e = aux_e;
	for (const auto &d : decomp) {
		switch (d.t) {
			case ML_UpT:
				e += ML_Up[d.br][d.used_mask][d.up][d.upr][d.i][d.j];
				break;
			case ML_BrT:
				e += ML_Br[d.br][d.used_mask][d.up][d.upr][d.i][d.j];
				break;
			case PT:
				if (parent == ET)
					e += SSScore(d.i, d.j);
				else
					e += MLSSScore(d.i, d.j);
				break;
			case CxFlT:
				e += CxFl[d.i][d.j];
				break;
			case CxMM3T:
				e += CxMM3[d.i][d.j];
				break;
			case CxMM5T:
				e += CxMM5[d.i][d.j];
				break;
			case ET:
				e += E[d.i];
				break;
		}
	}
	if (e < best) {
		best = e;
		best_decomp = decomp;
	}
}

void librnary::AsymmetryFolder::TraceE(std::stack<TState> &s) {
	// Manage current state.
	TState curr = s.top();
	int i = curr.i;
	s.pop();

	// Base case.
	if (i == 0)
		return;

	// Maintain best decomposition.
	vector<TState> best_decomp = {TState(i - 1)}; // Unpaired 3'.
	energy_t be = E[i - 1];
	for (int k = -1; k < i; ++k) {
		int decomp = k == -1 ? 0 : E[k];
		vector<TState> decomp_state;
		if (k != -1)
			decomp_state.push_back(TState(k));
		energy_t e = decomp + SSScore(k + 1, i);
		if (e < be) {
			be = e;
			best_decomp = {TState(k + 1, i)};
			best_decomp.insert(best_decomp.end(), decomp_state.begin(), decomp_state.end());
		}
		if (stacking) {
			if (k + 2 < i) {
				e = decomp + SSScore(k + 2, i) + em.FiveDangle(k + 2, i);
				if (e < be) {
					be = e;
					best_decomp = {TState(k + 2, i)};
					best_decomp.insert(best_decomp.end(), decomp_state.begin(), decomp_state.end());
				}
			}
			if (k + 1 < i - 1) {
				e = decomp + SSScore(k + 1, i - 1) + em.ThreeDangle(k + 1, i - 1);
				if (e < be) {
					be = e;
					best_decomp = {TState(k + 1, i - 1)};
					best_decomp.insert(best_decomp.end(), decomp_state.begin(), decomp_state.end());
				}
			}
			if (k + 2 < i - 1) {
				e = decomp + SSScore(k + 2, i - 1) + em.Mismatch(k + 2, i - 1);
				if (e < be) {
					be = e;
					best_decomp = {TState(k + 2, i - 1)};
					best_decomp.insert(best_decomp.end(), decomp_state.begin(), decomp_state.end());
				}
			}

			// Coaxial stack decompositions.
			for (int j = k + 1; j + 1 < i; ++j) {
				if (k + 1 < j) {
					e = decomp + em.FlushCoax(k + 1, j, j + 1, i)
						+ SSScore(k + 1, j) + SSScore(j + 1, i);
					if (e < be) {
						be = e;
						best_decomp = {TState(k + 1, j), TState(j + 1, i)};
						best_decomp.insert(best_decomp.end(), decomp_state.begin(), decomp_state.end());
					}
				}
				if (k + 2 < j - 1) {
					e = decomp + em.MismatchCoax(k + 2, j - 1, j + 1, i)
						+ SSScore(k + 2, j - 1) + SSScore(j + 1, i);
					if (e < be) {
						be = e;
						best_decomp = {TState(k + 2, j - 1), TState(j + 1, i)};
						best_decomp.insert(best_decomp.end(), decomp_state.begin(), decomp_state.end());
					}
				}
				if (j + 2 < i - 1 && k + 1 < j) {
					e = decomp + em.MismatchCoax(j + 2, i - 1, k + 1, j)
						+ SSScore(k + 1, j) + SSScore(j + 2, i - 1);
					if (e < be) {
						be = e;
						best_decomp = {TState(k + 1, j), TState(j + 2, i - 1)};
						best_decomp.insert(best_decomp.end(), decomp_state.begin(), decomp_state.end());
					}
				}
			}
		}

	}

	for (const TState &state : best_decomp)
		s.push(state);

}
void librnary::AsymmetryFolder::TraceP(std::stack<TState> &s) {
	TState curr = s.top();
	s.pop();
	int i = curr.i, j = curr.j;

	const int req_br = DefaultRequiredMultiInternalBranches();
	int up_lim = UnpairedGapLimit();

	// Hairpin. Also technically a base case.
	vector<TState> best_decomp;
	energy_t be = em.OneLoop(i, j);
	auto relax_decomp = [&](const vector<TState> &states, int aux_e) {
		return Relax(PT, be, best_decomp, states, aux_e);
	};
	// Multi-loops.
	for (int up = 0; up <= up_lim && i + 1 + up < j; ++up) {
		for (int upr = 0; upr <= up_lim && i + 1 + up < j - 1 - upr; ++upr) {
			// Base of a multi-loop cost.
			energy_t ml_base = em.Branch(i, j) + em.MLInit() + em.MLUnpairedCost() * up + em.MLBranchCost();
			energy_t ml_asym_base = ml_base + AsymScore(up, upr);
			relax_decomp({TState(ML_BrT, req_br - 1, 0, up, upr, i + 1 + up, j - 1)}, ml_asym_base);
			if (stacking) {
				// (._) left dangle.
				if (up >= 1) {
					relax_decomp({TState(ML_BrT, req_br - 1, 1, up, upr, i + 1 + up, j - 1)},
								 ml_asym_base + em.ClosingThreeDangle(i, j));
				}
				//(_.) right dangle.
				if (upr >= 1) {
					relax_decomp({TState(ML_BrT, req_br - 1, 2, up, upr, i + 1 + up, j - 1)},
								 ml_asym_base + em.ClosingFiveDangle(i, j));
				}
				//(._.) mismatch.
				if (upr >= 1 && up >= 1) {
					relax_decomp({TState(ML_BrT, req_br - 1, 3, up, upr, i + 1 + up, j - 1)},
								 ml_asym_base + em.ClosingMismatch(i, j));
				}
				// Coaxial stacks on the left.
				// The loop condition basically disallows ML_Br fragments without a possible branch.
				for (int k = i + 2; k + up + 1 < j - upr - 1; ++k) {
					// ((_)_)
					//    ^ <- k
					relax_decomp({TState(ML_BrT, req_br - 2, 0, up, upr, k + up + 1, j - 1), TState(i + 1, k)},
								 ml_base + AsymScore(0, up) + AsymScore(0, upr) + em.FlushCoax(i, j, i + 1, k));

					// (.(_)_.)
					if (upr >= 1 && i + 2 < k) {
						relax_decomp({TState(ML_BrT, req_br - 2, 2, up, upr, k + up + 1, j - 1), TState(i + 2, k)},
									 ml_base
										 + AsymScore(1, up) + AsymScore(1, upr) + em.MLUnpairedCost()
										 + em.MismatchCoax(i, j, i + 2, k));

					}
					// (.(_)._)
					if (up >= 1 && i + 2 < k) {
						relax_decomp({TState(ML_BrT, req_br - 2, 1, up, upr, k + up + 1, j - 1), TState(i + 2, k)},
									 ml_base
										 + AsymScore(1, up) + AsymScore(1, upr) + em.MLUnpairedCost()
										 + em.MismatchCoax(i + 2, k, i, j));
					}
				}
				// Coaxial stacks on the right.
				// Loop condition exists for the same reason as the previous coax loop.
				for (int k = j - 2; k - upr - 1 > i + 1 + up; --k) {
					// (_(_))
					//   ^ <- k
					relax_decomp({TState(ML_BrT, req_br - 2, 0, up, upr, i + 1 + up, k - 1), TState(k, j - 1)}, ml_base
						+ AsymScore(0, up) + AsymScore(0, upr) + em.FlushCoax(i, j, k, j - 1));

					// (._(_).)
					if (up >= 1 && k < j - 2) {
						relax_decomp({TState(ML_BrT, req_br - 2, 1, up, upr, i + 1 + up, k - 1), TState(k, j - 2)},
									 ml_base
										 + AsymScore(1, up) + AsymScore(1, upr) + em.MLUnpairedCost()
										 + em.MismatchCoax(i, j, k, j - 2));

					}
					// (_.(_).)
					if (upr >= 1 && k < j - 2) {
						relax_decomp({TState(ML_BrT, req_br - 2, 2, up, upr, i + 1 + up, k - 1), TState(k, j - 2)},
									 ml_base
										 + AsymScore(1, up) + AsymScore(1, upr) + em.MLUnpairedCost()
										 + em.MismatchCoax(k, j - 2, i, j));
					}
				}
			}

		}
	}
	// Be sure not to use relax_decomp here, because it will call MLSSScore.
	// Two loops for the two-loops (stack, bulge, or internal loop).
	for (int k = i + 1; k + 1 < j && (k - i - 1) <= max_twoloop_unpaired; ++k)
		for (int l = j - 1; l > k && (j - l - 1) + (k - i - 1) <= max_twoloop_unpaired; --l) {
			energy_t e = P[k][l] + em.TwoLoop(i, k, l, j);
			if (e < be) {
				be = e;
				best_decomp = {TState(k, l)};
			}
		}

	for (const TState &state : best_decomp)
		s.push(state);
}
void librnary::AsymmetryFolder::TraceCxFl(std::stack<TState> &s) {
	TState curr = s.top();
	s.pop();
	int i = curr.i, j = curr.j;

	vector<TState> best_decomp;
	energy_t be = em.MaxMFE();
	energy_t e;

	for (int k = i + 1; k + 1 < j; ++k) {
		e = em.FlushCoax(i, k, k + 1, j) + MLSSScore(i, k) + MLSSScore(k + 1, j);
		if (e < be) {
			be = e;
			best_decomp = {TState(i, k), TState(k + 1, j)};
		}
	}

	for (const TState &state : best_decomp)
		s.push(state);

}
void librnary::AsymmetryFolder::TraceCxMM3(std::stack<TState> &s) {
	TState curr = s.top();
	s.pop();
	int i = curr.i, j = curr.j;

	vector<TState> best_decomp;
	energy_t be = em.MaxMFE();
	auto relax_decomp = [&](const vector<TState> &states, int aux_e) {
		return Relax(CxMM3T, be, best_decomp, states, aux_e);
	};
	for (int k = i + 1; k + 2 < j; ++k) {
		relax_decomp({TState(i, k), TState(k + 2, j)}, em.MismatchCoax(k + 2, j, i, k));
	}

	for (const TState &state : best_decomp)
		s.push(state);
}
void librnary::AsymmetryFolder::TraceCxMM5(std::stack<TState> &s) {
	TState curr = s.top();
	s.pop();
	int i = curr.i, j = curr.j;

	vector<TState> best_decomp;
	energy_t be = em.MaxMFE();
	auto relax_decomp = [&](const vector<TState> &states, int aux_e) {
		return Relax(CxMM5T, be, best_decomp, states, aux_e);
	};
	for (int k = i + 1; k + 2 < j; ++k) {
		relax_decomp({TState(i, k), TState(k + 2, j)}, em.MismatchCoax(i, k, k + 2, j));
	}

	for (const TState &state : best_decomp)
		s.push(state);
}
void librnary::AsymmetryFolder::TraceMLUp(std::stack<TState> &s) {

	TState curr = s.top();
	s.pop();
	int i = curr.i, j = curr.j, br = curr.br, used_mask = curr.used_mask, up = curr.up, upr = curr.upr;
	int used5 = used_mask & 1;
//	int used3 = used_mask & 2;
	int br_prime = max(0, br - 1);

	vector<TState> best_decomp;
	energy_t be = em.MaxMFE();
	auto relax_decomp = [&](const vector<TState> &states, int aux_e) {
		return Relax(ML_UpT, be, best_decomp, states, aux_e);
	};

	if (br == 0 && (j - i + 1) == upr) {
		relax_decomp({}, upr * em.MLUnpairedCost() + AsymScore(up, upr));
	}
	// ...(_)xxx( -- must be followed by a branch
	// [i, k) covers no unpaired case naturally
	// Ensure that used5 is satisfied.
	for (int k = i + used5; (k - i) <= UnpairedGapLimit() && k < j; ++k) {
		relax_decomp({TState(ML_BrT, br_prime, used_mask, k - i, upr, k, j)},
					 em.MLUnpairedCost() * (k - i) + AsymScore(up, k - i));
	}


	for (const TState &state : best_decomp)
		s.push(state);

}
void librnary::AsymmetryFolder::TraceMLBr(std::stack<TState> &s) {
	TState curr = s.top();
	s.pop();
	int i = curr.i, j = curr.j, br = curr.br, used_mask = curr.used_mask, up = curr.up, upr = curr.upr;
	int used5 = used_mask & 1;
	int used3 = used_mask & 2;
	int br_prime = max(0, br - 1);

	vector<TState> best_decomp;
	energy_t be = em.MaxMFE();
	auto relax_decomp = [&](const vector<TState> &states, int aux_e) {
		return Relax(ML_BrT, be, best_decomp, states, aux_e);
	};

	for (int k = i + 1; k <= j; ++k) { // k <= j for end on feature cases.
		relax_decomp({TState(ML_UpT, br, used3, up, upr, k + 1, j), TState(i, k)}, 0);
		if (stacking) {
			// (_).{...}
			// Set min unpaired to 1 to force dangle.
			relax_decomp({TState(ML_UpT, br, used3 | 1, up, upr, k + 1, j), TState(i, k)}, em.ThreeDangle(i, k));
			if (up - used5 >= 1) {
				// .(_){...}
				relax_decomp({TState(ML_UpT, br, used3, up, upr, k + 1, j), TState(i, k)}, em.FiveDangle(i, k));
				// .(_).{...}
				// Set min unpaired to 1 to force dangle.
				relax_decomp({TState(ML_UpT, br, used3 | 1, up, upr, k + 1, j), TState(i, k)}, em.Mismatch(i, k));
			}
			// Coaxial stacks.
			// Carefully accounts for asymmetry between coax stacks, and br_prime
			// is replaced by 0.
			// (_)(_){...}
			relax_decomp({TState(CxFlT, i, k), TState(ML_UpT, br_prime, used3, 0, upr, k + 1, j)}, AsymScore(up, 0));

			// Needs to count the middle unpaired.
			// (_).(_).{...}
			// Makes sure to set min unpaired to 1, and need to include unpaired at k.
			relax_decomp({TState(CxMM3T, i, k), TState(ML_UpT, br_prime, used3 | 1, 1, upr, k + 1, j)},
						 AsymScore(up, 1) + em.MLUnpairedCost());
			// .(_).(_){...}
			if (up - used5 >= 1) {
				relax_decomp({TState(CxMM5T, i, k), TState(ML_UpT, br_prime, used3, 1, upr, k + 1, j)},
							 AsymScore(up, 1) + em.MLUnpairedCost());
			}
		}
	}

	for (const TState &state : best_decomp)
		s.push(state);
}

librnary::Matching librnary::AsymmetryFolder::Traceback() {
	const int N = static_cast<int>(rna.size());
	Matching m = EmptyMatching(static_cast<unsigned>(rna.size()));
	if (N == 0)
		return m;
	stack<TState> s;
	s.push(TState(N - 1));
	while (!s.empty()) {
		if (s.top().t == ET)
			TraceE(s);
		else if (s.top().t == PT) {
			m[s.top().i] = s.top().j;
			m[s.top().j] = s.top().i;
			TraceP(s);
		} else if (s.top().t == CxFlT)
			TraceCxFl(s);
		else if (s.top().t == CxMM3T)
			TraceCxMM3(s);
		else if (s.top().t == CxMM5T)
			TraceCxMM5(s);
		else if (s.top().t == ML_UpT)
			TraceMLUp(s);
		else // s.top().t == MLBrT
			TraceMLBr(s);
	}
	return m;
}


int librnary::AsymmetryFolder::Fold(const PrimeStructure &primary) {

	rna = primary;

	em.SetRNA(rna);

	const int N = static_cast<int>(rna.size());

	const int req_br = DefaultRequiredMultiInternalBranches();

	assert(req_br >= 2); // Otherwise we're not making multi-loops, and also seg-faults.

	// Ignore empty RNAs.
	if (N == 0)
		return 0;

	// Clear DP tables.
	E.clear();
	P.clear();
	if (stacking) {
		CxFl.clear();
		CxMM3.clear();
		CxMM5.clear();
	}

	// Resize DP tables.
	E.resize(rna.size(), 0);
	P.resize(rna.size(), VE(rna.size(), em.MaxMFE()));
	int up_lim = UnpairedGapLimit();
	unsigned up_sz = static_cast<unsigned>(up_lim + 1);
	ML_Up = Array2D<_4DVE>(static_cast<size_t>(req_br),
						   4,
						   _4DVE(up_sz + 1, VVVE(up_sz + 1, VVE(rna.size(), VE(rna.size(), em.MaxMFE())))));
	ML_Br = Array2D<_4DVE>(static_cast<size_t>(req_br),
						   4,
						   _4DVE(up_sz + 1, VVVE(up_sz + 1, VVE(rna.size(), VE(rna.size(), em.MaxMFE())))));
	if (stacking) {
		CxFl.resize(rna.size(), VE(rna.size(), em.MaxMFE()));
		CxMM5.resize(rna.size(), VE(rna.size(), em.MaxMFE()));
		CxMM3.resize(rna.size(), VE(rna.size(), em.MaxMFE()));
	}

	// Base case that allows the [i,i] fragment to end on a single unpaired.
	// Also allow [i,i-1] for empty fragments.
	for (int i = 1; i < N; ++i) { // We should never use fragments like [0,_] so start at i=1.
		// Number of forced unpaired.
		for (int up = 0; up <= up_lim; ++up) {
			// Either 5' unpaired is used, or closing 3' is unpaired is used, or neither are. Both can't be.
			ML_Up[0][0][up][1][i][i] = AsymScore(up, 1) + em.MLUnpairedCost();
			ML_Up[0][1][up][1][i][i] = AsymScore(up, 1) + em.MLUnpairedCost();
			ML_Up[0][2][up][1][i][i] = AsymScore(up, 1) + em.MLUnpairedCost();
			// Empty fragment case. No unpaired can be used.
			ML_Up[0][0][up][0][i][i - 1] = AsymScore(up, 0);
		}
	}

	for (int i = N - 1; i >= 0; --i) {
		for (int j = i + 1; j < N; ++j) {
			librnary::energy_t best;
			// Paired table.
			if (ValidPair(rna[i], rna[j]) &&
				(lonely_pairs || !MustBeLonelyPair(rna, i, j, em.MIN_HAIRPIN_UNPAIRED))) {
				best = em.OneLoop(i, j); // Hairpin.
				// Two loops for the two-loops (stack, bulge, or internal loosp).
				for (int k = i + 1; k + 1 < j && (k - i - 1) <= max_twoloop_unpaired; ++k)
					for (int l = j - 1; l > k && (j - l - 1) + (k - i - 1) <= max_twoloop_unpaired; --l)
						best = min(best, P[k][l] + em.TwoLoop(i, k, l, j));
				// Multi-loops.
				for (int up = 0; up <= up_lim && i + 1 + up < j; ++up) {
					for (int upr = 0; upr <= up_lim && i + 1 + up < j - 1 - upr; ++upr) {
						// Base of a multi-loop cost.
						energy_t ml_base = em.Branch(i, j) + em.MLInit() + em.MLUnpairedCost() * up + em.MLBranchCost();
						energy_t ml_asym_base = ml_base + AsymScore(up, upr);
						best = min(best, ML_Br[req_br - 1][0][up][upr][i + 1 + up][j - 1] + ml_asym_base);
						if (stacking) {
							// (._) left dangle.
							if (up >= 1) {
								best = min(best, ML_Br[req_br - 1][1][up][upr][i + 1 + up][j - 1]
									+ ml_asym_base + em.ClosingThreeDangle(i, j));
							}
							//(_.) right dangle.
							if (upr >= 1) {
								best = min(best, ML_Br[req_br - 1][2][up][upr][i + 1 + up][j - 1]
									+ ml_asym_base + em.ClosingFiveDangle(i, j));
							}
							//(._.) mismatch.
							if (upr >= 1 && up >= 1) {
								best = min(best, ML_Br[req_br - 1][3][up][upr][i + 1 + up][j - 1] + ml_asym_base
									+ em.ClosingMismatch(i, j));
							}
							// Coaxial stacks on the left.
							// The loop condition basically disallows ML_Br fragments without a possible branch.
							for (int k = i + 2; k + up + 1 < j - upr - 1; ++k) {
								// ((_)_)
								//    ^ <- k
								best = min(best, ML_Br[req_br - 2][0][up][upr][k + up + 1][j - 1] + ml_base
									+ AsymScore(0, up) + AsymScore(0, upr) + em.FlushCoax(i, j, i + 1, k)
									+ MLSSScore(i + 1, k));
								// (.(_)_.)
								if (upr >= 1 && i + 2 < k) {
									best = min(best, ML_Br[req_br - 2][2][up][upr][k + up + 1][j - 1] + ml_base
										+ AsymScore(1, up) + AsymScore(1, upr) + em.MLUnpairedCost()
										+ em.MismatchCoax(i, j, i + 2, k) + MLSSScore(i + 2, k));

								}
								// (.(_)._)
								if (up >= 1 && i + 2 < k) {
									best = min(best, ML_Br[req_br - 2][1][up][upr][k + up + 1][j - 1] + ml_base
										+ AsymScore(1, up) + AsymScore(1, upr) + em.MLUnpairedCost()
										+ em.MismatchCoax(i + 2, k, i, j) + MLSSScore(i + 2, k));
								}
							}
							// Coaxial stacks on the right.
							// Loop condition exists for the same reason as the previous coax loop.
							for (int k = j - 2; k - upr - 1 > i + 1 + up; --k) {
								// (_(_))
								//   ^ <- k
								best = min(best, ML_Br[req_br - 2][0][up][upr][i + 1 + up][k - 1] + ml_base
									+ AsymScore(0, up) + AsymScore(0, upr) + em.FlushCoax(i, j, k, j - 1)
									+ MLSSScore(k, j - 1));
								// (._(_).)
								if (up >= 1 && k < j - 2) {
									best = min(best, ML_Br[req_br - 2][1][up][upr][i + 1 + up][k - 1] + ml_base
										+ AsymScore(1, up) + AsymScore(1, upr) + em.MLUnpairedCost()
										+ em.MismatchCoax(i, j, k, j - 2) + MLSSScore(k, j - 2));

								}
								// (_.(_).)
								if (upr >= 1 && k < j - 2) {
									best = min(best, ML_Br[req_br - 2][2][up][upr][i + 1 + up][k - 1] + ml_base
										+ AsymScore(1, up) + AsymScore(1, upr) + em.MLUnpairedCost()
										+ em.MismatchCoax(k, j - 2, i, j) + MLSSScore(k, j - 2));
								}
							}
						}

					}
				}
				P[i][j] = best;
			}

			// Multi-loop fragment tables.
			// The fact that i>0 and j<N-1 are invalid multi-loop fragments is useful in other places.
			// It allows the tables to represent empty fragments using a range like this [i,i-1].
			if (i > 0 && j < N - 1) { // Avoid needless fragments that also cause invalid energy model calls.

				// Fill the coaxial stack tables if stacking is enabled.
				// These tables work like (_)(_) or (_).(_)
				//                        ^    ^    ^     ^
				//                        i    j    i     j
				// That is for CxMM5/3 the extra unpaired nucleotide is outside [i,j].
				for (int k = i + 1; k + 1 < j && stacking; ++k) {
					CxFl[i][j] = min(CxFl[i][j], em.FlushCoax(i, k, k + 1, j) + MLSSScore(i, k) + MLSSScore(k + 1, j));
					if (k + 2 < j) {
						CxMM5[i][j] = min(CxMM5[i][j], em.MismatchCoax(i, k, k + 2, j)
							+ MLSSScore(i, k) + MLSSScore(k + 2, j));
						CxMM3[i][j] = min(CxMM3[i][j], em.MismatchCoax(k + 2, j, i, k)
							+ MLSSScore(i, k) + MLSSScore(k + 2, j));
					}
				}

				for (int br = 0; br < req_br; ++br) {
					// The used bit-mask defines which nucleotides have been used.
					for (int used_mask = 0; used_mask < 4; ++used_mask) {
						for (int up = 0; up <= up_lim; ++up) {
							for (int upr = 0; upr <= up_lim; ++upr) {

								int br_prime = max(0, br - 1);

								// ML_Br.
								// The bitmask for just if the relevant unpaired nt in the 5' segment was used.
								// For ML_Br this is the 5' most unpaired in the preceding section of up nucs.
								// For ML_Up this is the 5' unpaired following the last branch.
								// Either 0 or 1.
								int used5 = used_mask & 1;
								// Bitmask for if closing branch used the most 3' unpaired nt.
								// Either 0 or 2.
								int used3 = used_mask & 2;
								best = em.MaxMFE();
								// Decompose into branch and then unpaired.
								// ...(_
								for (int k = i + 1; k <= j; ++k) { // k <= j for end on feature cases.
									best = min(best, MLSSScore(i, k) + ML_Up[br][used3][up][upr][k + 1][j]);
									if (stacking) {
										// (_).{...}
										// Set min unpaired to 1 to force dangle.
										best = min(best, em.ThreeDangle(i, k) + MLSSScore(i, k)
											+ ML_Up[br][used3 | 1][up][upr][k + 1][j]);
										if (up - used5 >= 1) {
											// .(_){...}
											best = min(best, em.FiveDangle(i, k) +
												MLSSScore(i, k) + ML_Up[br][used3][up][upr][k + 1][j]);
											// .(_).{...}
											// Set min unpaired to 1 to force dangle.
											best = min(best, em.Mismatch(i, k) +
												MLSSScore(i, k) + ML_Up[br][used3 | 1][up][upr][k + 1][j]);
										}
										// Coaxial stacks.
										// Carefully accounts for asymmetry between coax stacks, and br_prime
										// is replaced by 0.
										// (_)(_){...}
										best = min(best, CxFl[i][k] + ML_Up[br_prime][used3][0][upr][k + 1][j]
											+ AsymScore(up, 0));

										// Needs to count the middle unpaired.
										// (_).(_).{...}
										// Makes sure to set min unpaired to 1, and need to include unpaired at k.
										best = min(best, CxMM3[i][k] + ML_Up[br_prime][used3 | 1][1][upr][k + 1][j]
											+ AsymScore(up, 1) + em.MLUnpairedCost());
										// .(_).(_){...}
										if (up - used5 >= 1) {
											best = min(best, CxMM5[i][k] + ML_Up[br_prime][used3][1][upr][k + 1][j]
												+ AsymScore(up, 1) + em.MLUnpairedCost());
										}
									}
								}
								ML_Br[br][used_mask][up][upr][i][j] = best;

								// ML_Up stuff.
								best = em.MaxMFE();
								// All unpaired until the end of the fragment.
								// No need to check used_mask here since we always have >= 2 unpaired.
								if (br == 0 && (j - i + 1) == upr) {
									best = min(best, upr * em.MLUnpairedCost() + AsymScore(up, upr));
								}
								// ...(_)xxx( -- must be followed by a branch
								// [i, k) covers no unpaired case naturally
								// Ensure that used5 is satisfied.
								for (int k = i + used5; (k - i) <= up_lim && k < j; ++k) {
									best = min(best, em.MLUnpairedCost() * (k - i) + AsymScore(up, k - i)
										+ ML_Br[br_prime][used_mask][k - i][upr][k][j]);
								}
								ML_Up[br][used_mask][up][upr][i][j] = best;
							}
						}
					}
				}
			}
		}
	}

	for (int i = 1; i < N; ++i) {
		energy_t best = E[i - 1];
		for (int k = -1; k < i; ++k) {
			int decomp = k == -1 ? 0 : E[k];
			best = min(best, decomp + SSScore(k + 1, i));
			if (stacking) {
				if (k + 2 < i)
					best = min(best, decomp + SSScore(k + 2, i) + em.FiveDangle(k + 2, i));
				if (k + 1 < i - 1)
					best = min(best, decomp + SSScore(k + 1, i - 1) + em.ThreeDangle(k + 1, i - 1));
				if (k + 2 < i - 1)
					best = min(best, decomp + SSScore(k + 2, i - 1) + em.Mismatch(k + 2, i - 1));
				// Coaxial stack decompositions.
				for (int j = k + 1; j + 1 < i; ++j) {
					if (k + 1 < j)
						best = min(best, decomp + em.FlushCoax(k + 1, j, j + 1, i)
							+ SSScore(k + 1, j) + SSScore(j + 1, i));
					if (k + 2 < j - 1)
						best = min(best, decomp + em.MismatchCoax(k + 2, j - 1, j + 1, i)
							+ SSScore(k + 2, j - 1) + SSScore(j + 1, i));
					if (k + 1 < j && j + 2 < i - 1)
						best = min(best, decomp + em.MismatchCoax(j + 2, i - 1, k + 1, j)
							+ SSScore(k + 1, j) + SSScore(j + 2, i - 1));
				}
			}
		}
		E[i] = best;
	}

	return E[N - 1];
}
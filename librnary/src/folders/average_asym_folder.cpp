//
// Created by max on 8/11/16.
//

#include "folders/average_asym_folder.hpp"

#include <multi_loop.hpp>
#include <scorers/nn_scorer.hpp>
#include "models/nn_affine_model.hpp"
#include "scorers/average_asym_scorer.hpp"

using namespace std;

int librnary::AverageAsymmetryFolder::MaxMLBranches() const {
	return max_ml_branches;
}
void librnary::AverageAsymmetryFolder::SetMaxMLBranches(int v) {
	max_ml_branches = v;
}

int librnary::AverageAsymmetryFolder::MaxMLNonClosingAsym() const {
	return max_nonclosing_ml_sum_asym_;
}
void librnary::AverageAsymmetryFolder::SetMaxMLNonClosingAsym(int v) {
	max_nonclosing_ml_sum_asym_ = v;
}

int librnary::AverageAsymmetryFolder::BranchesUB() const {
	return min(static_cast<int>(rna.size()) / (2 + em.MIN_HAIRPIN_UNPAIRED) + 1, max_ml_branches);
}


int librnary::AverageAsymmetryFolder::NonClosingSumAsymmetryUB() const {
	// The max sum asymmetry comes from the fact that the worst multi-loop looks like this:
	// ((_)(_)_(_))
	// Hence the formula on the right.
	return min(max_nonclosing_ml_sum_asym_,
			   max(0, 2 * static_cast<int>(rna.size() - 2 - 3 * (2 + em.MIN_HAIRPIN_UNPAIRED))));
}

int librnary::AverageAsymmetryFolder::UnpairedGapUB() const {
	return min(static_cast<int>(rna.size()), max_unpaired_gap);
}


librnary::energy_t librnary::AverageAsymmetryFolder::SSScore(int i, int j) const {
	return em.Branch(i, j) + P[i][j];
}

librnary::energy_t librnary::AverageAsymmetryFolder::MLSSScore(int i, int j) const {
	return SSScore(i, j) + em.MLBranchCost();
}


int librnary::AverageAsymmetryFolder::UnpairedGap() const {
	return max_unpaired_gap;
}
void librnary::AverageAsymmetryFolder::SetUnpairedGap(unsigned v) {
	max_unpaired_gap = v;
}


void librnary::AverageAsymmetryFolder::SetMaxTwoLoop(unsigned max_up) {
	max_twoloop_unpaired = max_up;
}

int librnary::AverageAsymmetryFolder::MaxTwoLoop() const {
	return max_twoloop_unpaired;
}

void librnary::AverageAsymmetryFolder::SetLonelyPairs(bool v) {
	lonely_pairs = v;
}

bool librnary::AverageAsymmetryFolder::LonelyPairs() const {
	return lonely_pairs;
}

void librnary::AverageAsymmetryFolder::SetStacking(bool v) {
	stacking = v;
}

bool librnary::AverageAsymmetryFolder::Stacking() const {
	return stacking;
}


librnary::VVE librnary::AverageAsymmetryFolder::GetP() const {
	return P;
}

void librnary::AverageAsymmetryFolder::Relax(Table parent, int &best, vector<TState> &best_decomp,
											 const vector<TState> &decomp, int aux_e) const {
	int e = aux_e;
	for (const auto &d : decomp) {
		switch (d.t) {
			case ML_UpT:
				e += ML_Up[d.used_mask][d.br][d.sum_asym][d.up][d.upr][d.i][d.j];
				break;
			case ML_BrT:
				e += ML_Br[d.used_mask][d.br][d.sum_asym][d.up][d.upr][d.i][d.j];
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

void librnary::AverageAsymmetryFolder::TraceE(std::stack<TState> &s) {
	// Manage current state.
	TState curr = s.top();
	int i = curr.i;
	s.pop();

	// Base case.
	if (i == 0)
		return;

	// Maintain best decomposition.
	vector<TState> best_decomp = {TState(i - 1)}; // Unpaired 3'.
	int be = E[i - 1];
	for (int k = -1; k < i; ++k) {
		int decomp = k == -1 ? 0 : E[k];
		vector<TState> decomp_state;
		if (k != -1)
			decomp_state.push_back(TState(k));
		int e = decomp + SSScore(k + 1, i);
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

void librnary::AverageAsymmetryFolder::TraceP(std::stack<TState> &s) {
	TState curr = s.top();
	s.pop();
	int i = curr.i, j = curr.j;

	// Hairpin. Also technically a base case.
	vector<TState> best_decomp;
	int be = em.OneLoop(i, j);
	auto relax_decomp = [&](const vector<TState> &states, int aux_e) {
		return Relax(PT, be, best_decomp, states, aux_e);
	};

	int br_lim = BranchesUB();
	int sum_asym_lim = NonClosingSumAsymmetryUB();
	int up_lim = UnpairedGapUB();


	// Multi-loops with 3 branches for strain purposes.

	// (_(_)_(_)_)
	//   ^ ^ ^ ^
	//   a b c d
	AverageAsymmetryScorer scorer(em);
	scorer.SetRNA(rna);
	for (int a = i + 1; a - i - 1 <= up_lim && a < j; ++a) {
		for (int b = a + 1; b < j; ++b) {
			for (int c = b + 1; c - b - 1 <= up_lim && c < j; ++c) {
				for (int d = j - 1; d > c && j - d - 1 <= up_lim; --d) {
					int unpaired = (a - i - 1) + (c - b - 1) + (j - d - 1);
					int sum_asym = abs((a - i - 1) - (c - b - 1)) // Frist branch.
						+ abs((c - b - 1) - (j - d - 1)) // Second branch.
						+ abs((a - i - 1) - (j - d - 1)); // Closing branch.
					// Don't forget to remove the MLSSScore calls since relax will deal with them.
					int base_cost = em.MLInit() + em.MLClosureAsymCost(sum_asym, 3) + em.MLBranchCost()
						+ unpaired * em.MLUnpairedCost() + em.Branch(i, j);
					vector<TState> decomp = {TState(a, b), TState(c, d)};
					// No stacking interactions.
					relax_decomp(decomp, base_cost);

					// The stacking interactions.
					// Note: I am very proud, and also vaugly nauseated by, the following hack.
					if (stacking) {
						LoopRegion ml(i, j);
						ml.enclosed = {BondPair(a, b), BondPair(c, d)};
						relax_decomp(decomp, base_cost + scorer.OptimalMLStacking(ml));
					}
				}
			}
		}
	}

	// Multi-loops with more than 3 branches which could not be strained.
	for (int br = 3; br < br_lim; ++br) { // Note that this requires 3 or more for br.
		for (int sum_asym = 0; sum_asym <= sum_asym_lim; ++sum_asym) {
			for (int up = 0; up <= up_lim && i + 1 + up < j; ++up) {
				for (int upr = 0; upr <= up_lim && i + 1 + up < j - 1 - upr; ++upr) {
					// Base of a multi-loop cost.
					int ml_base = em.Branch(i, j) + em.MLInit()
						+ em.MLUnpairedCost() * up + em.MLBranchCost();
					relax_decomp({TState(ML_BrT, 0, br - 1, sum_asym, up, upr, i + 1 + up, j - 1)}, ml_base
						+ em.MLClosureAsymCost(sum_asym + abs(up - upr), br + 1));
					if (stacking) {
						// (._) left dangle.
						if (up >= 1) {
							relax_decomp({TState(ML_BrT, 1, br - 1, sum_asym, up, upr, i + 1 + up, j - 1)}, ml_base
								+ em.ClosingThreeDangle(i, j)
								+ em.MLClosureAsymCost(sum_asym + abs(up - upr), br + 1));
						}
						//(_.) right dangle.
						if (upr >= 1) {
							relax_decomp({TState(ML_BrT, 2, br - 1, sum_asym, up, upr, i + 1 + up, j - 1)},
										 ml_base + em.ClosingFiveDangle(i, j)
											 + em.MLClosureAsymCost(sum_asym + abs(up - upr), br + 1));
						}
						//(._.) mismatch.
						if (upr >= 1 && up >= 1) {
							relax_decomp({TState(ML_BrT, 3, br - 1, sum_asym, up, upr, i + 1 + up, j - 1)},
										 ml_base + em.MLClosureAsymCost(sum_asym + abs(up - upr), br + 1)
											 + em.ClosingMismatch(i, j));
						}
						// Coaxial stacks on the left.
						// The loop condition basically disallows ML_Br fragments without a possible branch.
						for (int k = i + 2; k + up + 1 < j - upr - 1; ++k) {
							// ((_)_)
							//    ^ <- k
							relax_decomp({TState(ML_BrT, 0, br - 2, sum_asym, up, upr, k + up + 1, j - 1),
										  TState(i + 1, k)}, ml_base
											 + em.MLClosureAsymCost(sum_asym + up + upr, br + 1)
											 + em.FlushCoax(i, j, i + 1, k));
							// (.(_)_.)
							if (upr >= 1 && i + 2 < k) {
								relax_decomp({TState(ML_BrT, 2, br - 2, sum_asym, up, upr, k + up + 1, j - 1),
											  TState(i + 2, k)},
											 ml_base + em.MLUnpairedCost()
												 + em.MLClosureAsymCost(sum_asym + abs(1 - up) + abs(1 - upr), br + 1)
												 + em.MismatchCoax(i, j, i + 2, k));

							}
							// (.(_)._)
							if (up >= 1 && i + 2 < k) {
								relax_decomp({TState(ML_BrT, 1, br - 2, sum_asym, up, upr, k + up + 1, j - 1),
											  TState(i + 2, k)},
											 ml_base + em.MLUnpairedCost()
												 + em.MLClosureAsymCost(sum_asym + abs(1 - up) + abs(1 - upr), br + 1)
												 + em.MismatchCoax(i + 2, k, i, j));
							}
						}
						// Coaxial stacks on the right.
						// Loop condition exists for the same reason as the previous coax loop.
						for (int k = j - 2; k - upr - 1 > i + 1 + up; --k) {
							// (_(_))
							//   ^ <- k
							relax_decomp({TState(ML_BrT, 0, br - 2, sum_asym, up, upr, i + 1 + up, k - 1),
										  TState(k, j - 1)}, ml_base
											 + em.MLClosureAsymCost(sum_asym + up + upr, br + 1)
											 + em.FlushCoax(i, j, k, j - 1));
							// (._(_).)
							if (up >= 1 && k < j - 2) {
								relax_decomp({TState(ML_BrT, 1, br - 2, sum_asym, up, upr, i + 1 + up, k - 1),
											  TState(k, j - 2)},
											 ml_base + em.MLUnpairedCost()
												 + em.MLClosureAsymCost(sum_asym + abs(1 - up) + abs(1 - upr), br + 1)
												 + em.MismatchCoax(i, j, k, j - 2));

							}
							// (_.(_).)
							if (upr >= 1 && k < j - 2) {
								relax_decomp({TState(ML_BrT, 2, br - 2, sum_asym, up, upr, i + 1 + up, k - 1),
											  TState(k, j - 2)},
											 ml_base + em.MLUnpairedCost()
												 + em.MLClosureAsymCost(sum_asym + abs(1 - up) + abs(1 - upr), br + 1)
												 + em.MismatchCoax(k, j - 2, i, j));
							}
						}
					}

				}
			}
		}
	}

	// Be sure not to use relax_decomp here, because it will call MLSSScore.
	// Two loops for the two-loops (stack, bulge, or internal loop).
	for (int k = i + 1; k + 1 < j && (k - i - 1) <= max_twoloop_unpaired; ++k)
		for (int l = j - 1; l > k && (j - l - 1) + (k - i - 1) <= max_twoloop_unpaired; --l) {
			int e = P[k][l] + em.TwoLoop(i, k, l, j);
			if (e < be) {
				be = e;
				best_decomp = {TState(k, l)};
			}
		}

	for (const TState &state : best_decomp)
		s.push(state);
}

void librnary::AverageAsymmetryFolder::TraceCxFl(std::stack<TState> &s) {
	TState curr = s.top();
	s.pop();
	int i = curr.i, j = curr.j;

	vector<TState> best_decomp;
	int be = em.MaxMFE();
	int e;

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

void librnary::AverageAsymmetryFolder::TraceCxMM3(std::stack<TState> &s) {
	TState curr = s.top();
	s.pop();
	int i = curr.i, j = curr.j;

	vector<TState> best_decomp;
	int be = em.MaxMFE();
	auto relax_decomp = [&](const vector<TState> &states, int aux_e) {
		return Relax(CxMM3T, be, best_decomp, states, aux_e);
	};
	for (int k = i + 1; k + 2 < j; ++k) {
		relax_decomp({TState(i, k), TState(k + 2, j)}, em.MismatchCoax(k + 2, j, i, k));
	}

	for (const TState &state : best_decomp)
		s.push(state);
}

void librnary::AverageAsymmetryFolder::TraceCxMM5(std::stack<TState> &s) {
	TState curr = s.top();
	s.pop();
	int i = curr.i, j = curr.j;

	vector<TState> best_decomp;
	int be = em.MaxMFE();
	auto relax_decomp = [&](const vector<TState> &states, int aux_e) {
		return Relax(CxMM5T, be, best_decomp, states, aux_e);
	};
	for (int k = i + 1; k + 2 < j; ++k) {
		relax_decomp({TState(i, k), TState(k + 2, j)}, em.MismatchCoax(i, k, k + 2, j));
	}

	for (const TState &state : best_decomp)
		s.push(state);
}

void librnary::AverageAsymmetryFolder::TraceMLUp(std::stack<TState> &s) {

	TState curr = s.top();
	s.pop();
	int i = curr.i, j = curr.j, br = curr.br, sum_asym = curr.sum_asym,
		used_mask = curr.used_mask, up = curr.up, upr = curr.upr;
	int used5 = used_mask & 1;

	int up_lim = UnpairedGapUB();

	vector<TState> best_decomp;
	int be = em.MaxMFE();
	auto relax_decomp = [&](const vector<TState> &states, int aux_e) {
		return Relax(ML_UpT, be, best_decomp, states, aux_e);
	};


	// All unpaired until the end of the fragment.
	// No need to check used_mask here since we always have >= 2 unpaired.
	if (br == 0 && (j - i + 1) == upr && sum_asym == abs(up - upr)) {
		relax_decomp({}, upr * em.MLUnpairedCost());
	}

	if (br >= 1) {
		// ...(_)xxx( -- must be followed by a branch
		// [i, k) covers no unpaired case naturally
		// Ensure that used5 is satisfied.
		for (int k = i + used5; (k - i) <= up_lim && k < j; ++k) {
			if (sum_asym >= abs(up - (k - i)))
				relax_decomp({TState(ML_BrT, used_mask, br - 1, sum_asym - abs(up - (k - i)), k - i, upr, k, j)},
							 em.MLUnpairedCost() * (k - i));
		}
	}


	for (const TState &state : best_decomp)
		s.push(state);

}
void librnary::AverageAsymmetryFolder::TraceMLBr(std::stack<TState> &s) {
	TState curr = s.top();
	s.pop();
	int i = curr.i, j = curr.j, br = curr.br, used_mask = curr.used_mask, up = curr.up,
		upr = curr.upr, sum_asym = curr.sum_asym;
	int used5 = used_mask & 1;
	int used3 = used_mask & 2;

	vector<TState> best_decomp;
	int be = em.MaxMFE();
	auto relax_decomp = [&](const vector<TState> &states, int aux_e) {
		return Relax(ML_BrT, be, best_decomp, states, aux_e);
	};

	// Decompose into branch and then unpaired.
	// ...(_
	for (int k = i + 1; k <= j; ++k) { // k <= j for end on feature cases.
		relax_decomp({TState(ML_UpT, used3, br, sum_asym, up, upr, k + 1, j), TState(i, k)}, 0);
		if (stacking) {
			// (_).{...}
			// Set min unpaired to 1 to force dangle.
			relax_decomp({TState(ML_UpT, used3 | 1, br, sum_asym, up, upr, k + 1, j), TState(i, k)},
						 em.ThreeDangle(i, k));
			if (up - used5 >= 1) {
				// .(_){...}
				relax_decomp({TState(ML_UpT, used3, br, sum_asym, up, upr, k + 1, j), TState(i, k)},
							 em.FiveDangle(i, k));
				// .(_).{...}
				// Set min unpaired to 1 to force dangle.
				relax_decomp({TState(ML_UpT, used3 | 1, br, sum_asym, up, upr, k + 1, j), TState(i, k)},
							 em.Mismatch(i, k));
			}
			// Coaxial stacks.
			if (br >= 1) {
				// Carefully accounts for asymmetry between coax stacks,
				// (_)(_){...}
				if (sum_asym >= up)
					relax_decomp({TState(ML_UpT, used3, br - 1, sum_asym - up, 0, upr, k + 1, j), TState(CxFlT, i, k)},
								 0);

				// Needs to count the middle unpaired.
				// (_).(_).{...}
				// Makes sure to set min unpaired to 1, and need
				// to include unpaired at k.
				if (sum_asym >= abs(up - 1)) {
					relax_decomp({TState(ML_UpT, used3 | 1, br - 1, sum_asym - abs(up - 1), 1, upr, k + 1, j),
								  TState(CxMM3T, i, k)}, em.MLUnpairedCost());
					// .(_).(_){...}
					if (up - used5 >= 1) {
						relax_decomp({TState(ML_UpT, used3, br - 1, sum_asym - abs(up - 1), 1, upr, k + 1, j),
									  TState(CxMM5T, i, k)}, em.MLUnpairedCost());
					}
				}
			}
		}
	}

	for (const TState &state : best_decomp)
		s.push(state);
}

librnary::Matching librnary::AverageAsymmetryFolder::Traceback() {
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

int librnary::AverageAsymmetryFolder::Fold(const PrimeStructure &primary) {

	rna = primary;

	em.SetRNA(rna);

	const int N = static_cast<int>(rna.size());

	if (N == 0)
		return 0;

	// Clear DP tables.
	E.clear();
	P.clear();
	ML_Up.clear();
	ML_Br.clear();
	CxFl.clear();
	CxMM3.clear();
	CxMM5.clear();

	// Resize DP tables.
	E.resize(rna.size(), 0);
	P.resize(rna.size(), VI(rna.size(), em.MaxMFE()));
	int up_lim = UnpairedGapUB();
	// The maximum number of multi-loop branches we need to consider.
	int br_lim = BranchesUB();
	// Avoid the case where no multi-loops are allowed to avoid crashes or dealing with many special cases.
	if (br_lim <= 1) {
		br_lim = 2;
	}
	int sum_asym_lim = NonClosingSumAsymmetryUB();
	size_t up_sz = static_cast<size_t>(up_lim);
	size_t max_br = static_cast<size_t>(br_lim);
	size_t max_sum_asym = static_cast<size_t>(sum_asym_lim);
	// Note, we can only store up to max_br-1 because the closing and first branch is always done in the P table.
	ML_Up.resize(4, _6DVE(max_br - 1, _5DVE(max_sum_asym + 1, _4DVE(up_sz + 1,
																	VVVE(up_sz + 1,
																		 VVE(rna.size(),
																			 VE(rna.size(), em.MaxMFE())))))));
	ML_Br.resize(4, _6DVE(max_br - 1, _5DVE(max_sum_asym + 1, _4DVE(up_sz + 1,
																	VVVE(up_sz + 1,
																		 VVE(rna.size(),
																			 VE(rna.size(), em.MaxMFE())))))));
	if (stacking) {
		CxFl.resize(rna.size(), VI(rna.size(), em.MaxMFE()));
		CxMM5.resize(rna.size(), VI(rna.size(), em.MaxMFE()));
		CxMM3.resize(rna.size(), VI(rna.size(), em.MaxMFE()));
	}

	// Base case that allows the [i,i] fragment to end on a single unpaired.
	// Also allow [i,i-1] for empty fragments.
	for (int i = 1; i < N; ++i) { // We should never use fragments like [0,_] so start at i=1.
		// Number of forced unpaired.
		for (int up = 0; up <= up_lim; ++up) {
			// Either 5' unpaired is used, or closing 3' is unpaired is used, or neither are. Both can't be.
			// Also ensure the cases have the right amount of sum asymmetry.
			for (int bs = 0; bs < 3; ++bs) {
				if (up_lim >= 1 && sum_asym_lim >= abs(up - 1))
					ML_Up[bs][0][abs(up - 1)][up][1][i][i] = em.MLUnpairedCost();
			}
			// Empty fragment case. No unpaired can be used. Ensure min sum asym = 0.
			if (sum_asym_lim >= up)
				ML_Up[0][0][up][up][0][i][i - 1] = 0;
		}
	}

	for (int i = N - 1; i >= 0; --i) {
		for (int j = i + 1; j < N; ++j) {
			energy_t best;
			// Paired table.
			if (ValidPair(rna[i], rna[j]) &&
				(lonely_pairs || !MustBeLonelyPair(rna, i, j, em.MIN_HAIRPIN_UNPAIRED))) {
				best = em.OneLoop(i, j); // Hairpin.
				// Two loops for the two-loops (stack, bulge, or internal loosp).
				for (int k = i + 1; k + 1 < j && (k - i - 1) <= max_twoloop_unpaired; ++k)
					for (int l = j - 1; l > k && (j - l - 1) + (k - i - 1) <= max_twoloop_unpaired; --l)
						best = min(best, P[k][l] + em.TwoLoop(i, k, l, j));
				// Multi-loops.

				// Multi-loops with 3 branches for strain purposes.

				// (_(_)_(_)_)
				//   ^ ^ ^ ^
				//   a b c d
				if (br_lim >= 3) {
					AverageAsymmetryScorer scorer(em);
					scorer.SetRNA(rna);
					for (int a = i + 1; a - i - 1 <= up_lim && a < j; ++a) {
						for (int b = a + 1; b < j; ++b) {
							for (int c = b + 1; c - b - 1 <= up_lim && c < j; ++c) {
								for (int d = j - 1; d > c && j - d - 1 <= up_lim; --d) {
									int unpaired = (a - i - 1) + (c - b - 1) + (j - d - 1);
									int sum_asym = abs((a - i - 1) - (c - b - 1)) // Frist branch.
										+ abs((c - b - 1) - (j - d - 1)) // Second branch.
										+ abs((a - i - 1) - (j - d - 1)); // Closing branch.
									int base_cost = em.MLInit()
										+ em.MLClosureAsymCost(sum_asym, 3) + unpaired * em.MLUnpairedCost()
										+ em.Branch(i, j) + MLSSScore(a, b) + MLSSScore(c, d) + em.MLBranchCost();
									// No stacking interactions.
									best = min(best, base_cost);

									// The stacking interactions.
									// Note: I am very proud, and also vaugly nauseated by, the following hack.
									if (stacking) {
										LoopRegion ml(i, j);
										ml.enclosed = {BondPair(a, b), BondPair(c, d)};
										best = min(best, base_cost + scorer.OptimalMLStacking(ml));
									}
								}
							}
						}
					}
				}

				// Multi-loops with more than 3 internal branches which could not be strained.
				for (int br = 3; br < br_lim; ++br) { // Note that br >= 3.
					for (int sum_asym = 0; sum_asym <= sum_asym_lim; ++sum_asym) {
						for (int up = 0; up <= up_lim && i + 1 + up < j; ++up) {
							for (int upr = 0; upr <= up_lim && i + 1 + up < j - 1 - upr; ++upr) {
								// Base of a multi-loop cost.
								int ml_base = em.Branch(i, j) + em.MLInit()
									+ em.MLUnpairedCost() * up + em.MLBranchCost();
								best = min(best, ML_Br[0][br - 1][sum_asym][up][upr][i + 1 + up][j - 1] + ml_base
									+ em.MLClosureAsymCost(sum_asym + abs(up - upr), br + 1));
								if (stacking) {
									// (._) left dangle.
									if (up >= 1) {
										best = min(best, ML_Br[1][br - 1][sum_asym][up][upr][i + 1 + up][j - 1]
											+ ml_base + em.ClosingThreeDangle(i, j)
											+ em.MLClosureAsymCost(sum_asym + abs(up - upr), br + 1));
									}
									//(_.) right dangle.
									if (upr >= 1) {
										best = min(best, ML_Br[2][br - 1][sum_asym][up][upr][i + 1 + up][j - 1]
											+ ml_base + em.ClosingFiveDangle(i, j)
											+ em.MLClosureAsymCost(sum_asym + abs(up - upr), br + 1));
									}
									//(._.) mismatch.
									if (upr >= 1 && up >= 1) {
										best = min(best, ML_Br[3][br - 1][sum_asym][up][upr][i + 1 + up][j - 1]
											+ ml_base + em.MLClosureAsymCost(sum_asym + abs(up - upr), br + 1)
											+ em.ClosingMismatch(i, j));
									}
									// Coaxial stacks on the left.
									// The loop condition basically disallows ML_Br fragments without a possible branch.
									for (int k = i + 2; k + up + 1 < j - upr - 1; ++k) {
										// ((_)_)
										//    ^ <- k
										best =
											min(best, ML_Br[0][br - 2][sum_asym][up][upr][k + up + 1][j - 1] + ml_base
												+ em.MLClosureAsymCost(sum_asym + up + upr, br + 1)
												+ em.FlushCoax(i, j, i + 1, k)
												+ MLSSScore(i + 1, k));
										// (.(_)_.)
										if (upr >= 1 && i + 2 < k) {
											best = min(best, ML_Br[2][br - 2][sum_asym][up][upr][k + up + 1][j - 1]
												+ ml_base + em.MLUnpairedCost()
												+ em.MLClosureAsymCost(sum_asym + abs(1 - up) + abs(1 - upr), br + 1)
												+ em.MismatchCoax(i, j, i + 2, k) + MLSSScore(i + 2, k));

										}
										// (.(_)._)
										if (up >= 1 && i + 2 < k) {
											best = min(best, ML_Br[1][br - 2][sum_asym][up][upr][k + up + 1][j - 1]
												+ ml_base + em.MLUnpairedCost()
												+ em.MLClosureAsymCost(sum_asym + abs(1 - up) + abs(1 - upr), br + 1)
												+ em.MismatchCoax(i + 2, k, i, j) + MLSSScore(i + 2, k));
										}
									}
									// Coaxial stacks on the right.
									// Loop condition exists for the same reason as the previous coax loop.
									for (int k = j - 2; k - upr - 1 > i + 1 + up; --k) {
										// (_(_))
										//   ^ <- k
										best =
											min(best, ML_Br[0][br - 2][sum_asym][up][upr][i + 1 + up][k - 1] + ml_base
												+ em.MLClosureAsymCost(sum_asym + up + upr, br + 1)
												+ em.FlushCoax(i, j, k, j - 1) + MLSSScore(k, j - 1));
										// (._(_).)
										if (up >= 1 && k < j - 2) {
											best = min(best, ML_Br[1][br - 2][sum_asym][up][upr][i + 1 + up][k - 1]
												+ ml_base + em.MLUnpairedCost()
												+ em.MLClosureAsymCost(sum_asym + abs(1 - up) + abs(1 - upr), br + 1)
												+ em.MismatchCoax(i, j, k, j - 2) + MLSSScore(k, j - 2));

										}
										// (_.(_).)
										if (upr >= 1 && k < j - 2) {
											best = min(best, ML_Br[2][br - 2][sum_asym][up][upr][i + 1 + up][k - 1]
												+ ml_base + em.MLUnpairedCost()
												+ em.MLClosureAsymCost(sum_asym + abs(1 - up) + abs(1 - upr), br + 1)
												+ em.MismatchCoax(k, j - 2, i, j) + MLSSScore(k, j - 2));
										}
									}
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


				// The used bit-mask defines which nucleotides have been used.
				for (int used_mask = 0; used_mask < 4; ++used_mask) {
					// Because (internal branches <= multi-loop branches - 1) and P places a first branch,
					// this loop only needs to go < br_lim-1.
					for (int br = 0; br < br_lim - 1; ++br) { // See P DP fill for explanation of -1.
						for (int sum_asym = 0; sum_asym <= sum_asym_lim; ++sum_asym) {
							for (int up = 0; up <= up_lim; ++up) {
								for (int upr = 0; upr <= up_lim; ++upr) {
									// The bitmask for just if the relevant unpaired nt in the 5' segment was used.
									// For ML_Br this is the 5' most unpaired in the preceding section of up nucs.
									// For ML_Up this is the 5' unpaired following the last branch.
									// Either 0 or 1.
									int used5 = used_mask & 1;
									// Bitmask for if closing branch used the most 3' unpaired nt.
									// Either 0 or 2.
									int used3 = used_mask & 2;

									// ML_Br.
									// Don't allow fragments that start with a branch, but contain no branches.

									best = em.MaxMFE();
									// Decompose into branch and then unpaired.
									// ...(_
									for (int k = i + 1; k <= j; ++k) { // k <= j for end on feature cases.
										best = min(best, MLSSScore(i, k)
											+ ML_Up[used3][br][sum_asym][up][upr][k + 1][j]);
										if (stacking) {
											// (_).{...}
											// Set min unpaired to 1 to force dangle.
											best = min(best, em.ThreeDangle(i, k) + MLSSScore(i, k)
												+ ML_Up[used3 | 1][br][sum_asym][up][upr][k + 1][j]);
											if (up - used5 >= 1) {
												// .(_){...}
												best = min(best, em.FiveDangle(i, k) + MLSSScore(i, k)
													+ ML_Up[used3][br][sum_asym][up][upr][k + 1][j]);
												// .(_).{...}
												// Set min unpaired to 1 to force dangle.
												best = min(best, em.Mismatch(i, k) + MLSSScore(i, k)
													+ ML_Up[used3 | 1][br][sum_asym][up][upr][k + 1][j]);
											}
											// Coaxial stacks.
											if (br >= 1) {
												// Carefully accounts for asymmetry between coax stacks,
												// and br_prime is replaced by 0.
												// (_)(_){...}
												if (sum_asym >= up)
													best = min(best, CxFl[i][k]
														+ ML_Up[used3][br - 1][sum_asym - up][0][upr][k + 1][j]);

												// Needs to count the middle unpaired.
												// (_).(_).{...}
												// Makes sure to set min unpaired to 1, and need
												// to include unpaired at k.
												if (sum_asym >= abs(up - 1)) {
													best = min(best, CxMM3[i][k]
														+ ML_Up[used3 | 1][br - 1][sum_asym - abs(up - 1)][1][upr][k
															+ 1][j] + em.MLUnpairedCost());
													// .(_).(_){...}
													if (up - used5 >= 1) {
														best =
															min(best, CxMM5[i][k] + ML_Up[used3][br - 1][sum_asym -
																abs(up - 1)][1][upr][k + 1][j]
																+ em.MLUnpairedCost());
													}
												}
											}
										}
									}

									ML_Br[used_mask][br][sum_asym][up][upr][i][j] = best;


									// ML_Up stuff.
									best = em.MaxMFE();
									// All unpaired until the end of the fragment.
									// No need to check used_mask here since we always have >= 2 unpaired.
									if (br == 0 && (j - i + 1) == upr && sum_asym == abs(up - upr)) {
										best = min(best, upr * em.MLUnpairedCost());
									}

									if (br >= 1) {
										// ...(_)xxx( -- must be followed by a branch
										// [i, k) covers no unpaired case naturally
										// Ensure that used5 is satisfied.
										for (int k = i + used5; (k - i) <= up_lim && k < j; ++k) {
											if (sum_asym >= abs(up - (k - i)))
												best = min(best, em.MLUnpairedCost() * (k - i)
													+ ML_Br[used_mask][br - 1][sum_asym - abs(up - (k - i))][k
														- i][upr][k][j]);
										}
									}
									ML_Up[used_mask][br][sum_asym][up][upr][i][j] = best;
								}
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
			energy_t decomp = k == -1 ? 0 : E[k];
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
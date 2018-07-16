//
// Created by max on 6/22/16.
//

#include "folders/upbr_folder.hpp"

using namespace std;

int librnary::UpBrFolder::SSScore(int i, int j) {
	return em.Branch(i, j) + P[i][j];
}

void librnary::UpBrFolder::TraceE(stack<TState> &s) {
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

	for (const TState &state : best_decomp)
		s.push(state);
}

void librnary::UpBrFolder::TraceP(stack<TState> &s) {
	TState curr = s.top();
	s.pop();
	int i = curr.i, j = curr.j;

	// Hairpin. Also technically a base case.
	vector<TState> best_decomp; // Unpaired 3'.
	int be = em.OneLoop(i, j);
	int e;
	for (int up = 0; up < min(max_multi_unpaired + 1, j - i); ++up) { // Multiloops.
		for (int br = 2; br <= min(max_multi_branches, (j - i - 1 - up) / 2); ++br) {
			int init = em.MLClosure(up, br + 1) + em.Branch(i, j);
			e = ML[br][up][i + 1][j - 1] + init;
			if (e < be) {
				be = e;
				best_decomp = {TState(br, up, i + 1, j - 1)};
			}
			if (i + 2 < j - 1 && up >= 1) { // Left dangle.
				e = ML[br][up - 1][i + 2][j - 1] + init + em.ClosingThreeDangle(i, j);
				if (e < be) {
					be = e;
					best_decomp = {TState(br, up - 1, i + 2, j - 1)};
				}
			}
			if (i + 1 < j - 2 && up >= 1) { // Right dangle.
				e = ML[br][up - 1][i + 1][j - 2] + init + em.ClosingFiveDangle(i, j);
				if (e < be) {
					be = e;
					best_decomp = {TState(br, up - 1, i + 1, j - 2)};
				}
			}
			if (i + 2 < j - 2 && up >= 2) { // Mismatch.
				e = ML[br][up - 2][i + 2][j - 2] + init + em.ClosingMismatch(i, j);
				if (e < be) {
					be = e;
					best_decomp = {TState(br, up - 2, i + 2, j - 2)};
				}
			}


			// Coaxial stack.
			for (int k = i + 1; k < j; ++k) {
				// Five prime.
				// ((_)_)
				//    ^ <- k
				if (k + 1 < j - 1 && i + 1 < k) {
					e = ML[br - 1][up][k + 1][j - 1] + init
						+ em.FlushCoax(i, j, i + 1, k) +
						SSScore(i + 1, k);
					if (e < be) {
						be = e;
						best_decomp = {TState(br - 1, up, k + 1, j - 1), TState(i + 1, k)};
					}
				}
				if (up >= 2) {
					// (.(_)_.)
					if (i + 2 < k && k + 1 < j - 2) {
						e = ML[br - 1][up - 2][k + 1][j - 2] + init
							+ em.MismatchCoax(i, j, i + 2, k) +
							SSScore(i + 2, k);
						if (e < be) {
							be = e;
							best_decomp = {TState(br - 1, up - 2, k + 1, j - 2), TState(i + 2, k)};
						}
					}
					// (.(_)._)
					if (i + 2 < k && k + 2 < j) {
						e = ML[br - 1][up - 2][k + 2][j - 1] + init
							+ em.MismatchCoax(i + 2, k, i, j) +
							SSScore(i + 2, k);
						if (e < be) {
							be = e;
							best_decomp = {TState(br - 1, up - 2, k + 2, j - 1), TState(i + 2, k)};
						}
					}
				}

				// Three prime.
				// (_(_))
				//   ^ <- k
				if (i + 1 < k - 1 && k < j - 1) {
					e = ML[br - 1][up][i + 1][k - 1] + init
						+ em.FlushCoax(i, j, k, j - 1) +
						SSScore(k, j - 1);
					if (e < be) {
						be = e;
						best_decomp = {TState(br - 1, up, i + 1, k - 1), TState(k, j - 1)};
					}
				}
				if (up >= 2) {
					// (._(_).)
					if (k < j - 2 && i + 2 < k - 1) {
						e = ML[br - 1][up - 2][i + 2][k - 1] + init
							+ em.MismatchCoax(i, j, k, j - 2) +
							SSScore(k, j - 2);
						if (e < be) {
							be = e;
							best_decomp = {TState(br - 1, up - 2, i + 2, k - 1), TState(k, j - 2)};
						}
					}
					// (_.(_).)
					if (k < j - 2 && i + 1 < k - 2) {
						e = ML[br - 1][up - 2][i + 1][k - 2] + init
							+ em.MismatchCoax(k, j - 2, i, j) +
							SSScore(k, j - 2);
						if (e < be) {
							be = e;
							best_decomp = {TState(br - 1, up - 2, i + 1, k - 2), TState(k, j - 2)};
						}
					}
				}
			}
		}
	}
	// Two loops for the two-loops (stack, bulge, or internal loop).
	for (int k = i + 1; k + 1 < j && (k - i - 1) <= max_twoloop_unpaired; ++k)
		for (int l = j - 1; l > k && (j - l - 1) + (k - i - 1) <= max_twoloop_unpaired; --l) {
			e = P[k][l] + em.TwoLoop(i, k, l, j);
			if (e < be) {
				be = e;
				best_decomp = {TState(k, l)};
			}
		}

	for (const TState &state : best_decomp)
		s.push(state);
}

void librnary::UpBrFolder::TraceCxFl(stack<TState> &s) {
	TState curr = s.top();
	s.pop();
	int i = curr.i, j = curr.j;

	vector<TState> best_decomp;
	int be = em.MaxMFE();
	int e;

	for (int k = i + 1; k + 1 < j; ++k) {
		e = em.FlushCoax(i, k, k + 1, j) + SSScore(i, k) + SSScore(k + 1, j);
		if (e < be) {
			be = e;
			best_decomp = {TState(i, k), TState(k + 1, j)};
		}
	}

	for (const TState &state : best_decomp)
		s.push(state);

}

void librnary::UpBrFolder::TraceCxMM(stack<TState> &s) {
	TState curr = s.top();
	s.pop();
	int i = curr.i, j = curr.j;

	vector<TState> best_decomp;
	int be = em.MaxMFE();
	int e;

	for (int k = i + 1; k + 1 < j; ++k) {
		if (i + 1 < k - 1) {
			e = em.MismatchCoax(i + 1, k - 1, k + 1, j) + SSScore(i + 1, k - 1) +
				SSScore(k + 1, j);
			if (e < be) {
				be = e;
				best_decomp = {TState(i + 1, k - 1), TState(k + 1, j)};
			}
		}
		if (k + 2 < j - 1) {
			e = em.MismatchCoax(k + 2, j - 1, i, k) + SSScore(i, k) + SSScore(k + 2, j - 1);
			if (e < be) {
				be = e;
				best_decomp = {TState(i, k), TState(k + 2, j - 1)};
			}
		}
	}

	for (const TState &state : best_decomp)
		s.push(state);

}

void librnary::UpBrFolder::TraceML(stack<TState> &s) {
	TState curr = s.top();
	s.pop();
	int br = curr.b, up = curr.up, i = curr.i, j = curr.j;

	// Base case.
	if (i == j)
		return;

	vector<TState> best_decomp;
	int be = em.MaxMFE();

	if (up > 0) { // Single stranded.
		best_decomp = {TState(br, up - 1, i, j - 1)};
		be = ML[br][up - 1][i][j - 1];
	}

	int e;
	if (br == 1) { // End on branch cases.
		if (up == 0) {
			e = SSScore(i, j);
			if (e < be) {
				be = e;
				best_decomp = {TState(i, j)};
			}
		}
		if (i + 1 < j && up == 1) {
			e = SSScore(i + 1, j) + em.FiveDangle(i + 1, j);
			if (e < be) {
				be = e;
				best_decomp = {TState(i + 1, j)};
			}
			e = SSScore(i, j - 1) + em.ThreeDangle(i, j - 1);
			if (e < be) {
				be = e;
				best_decomp = {TState(i, j - 1)};
			}
		}
		if (i + 1 < j - 1 && up == 2) {
			e = SSScore(i + 1, j - 1) + em.Mismatch(i + 1, j - 1);
			if (e < be) {
				be = e;
				best_decomp = {TState(i + 1, j - 1)};
			}
		}
	}

	// End on coaxial stack.
	if (br == 2) {
		if (up == 2) {
			e = CxMM[i][j];
			if (e < be) {
				be = e;
				best_decomp = {TState(CxMMT, i, j)};
			}
		} else if (up == 0) {
			e = CxFl[i][j];
			if (e < be) {
				be = e;
				best_decomp = {TState(CxFlT, i, j)};
			}
		}
	}

	for (int k = i; k + 2 <= j; ++k) { // Try all decompositions into 5' ML fragment and 3' branch.
		if (br >= 1) {
			e = ML[br - 1][up][i][k] + SSScore(k + 1, j);
			if (e < be) {
				be = e;
				best_decomp = {TState(br - 1, up, i, k), TState(k + 1, j)};
			}
			if (k + 2 < j && up >= 1) {
				e = ML[br - 1][up - 1][i][k] + SSScore(k + 2, j) + em.FiveDangle(k + 2, j);
				if (e < be) {
					be = e;
					best_decomp = {TState(br - 1, up - 1, i, k), TState(k + 2, j)};
				}
			}
			if (k + 1 < j - 1 && up >= 1) {
				e = ML[br - 1][up - 1][i][k] + SSScore(k + 1, j - 1) +
					em.ThreeDangle(k + 1, j - 1);
				if (e < be) {
					be = e;
					best_decomp = {TState(br - 1, up - 1, i, k), TState(k + 1, j - 1)};
				}
			}
			if (k + 2 < j - 1 && up >= 2) {
				e = ML[br - 1][up - 2][i][k] + SSScore(k + 2, j - 1) +
					em.Mismatch(k + 2, j - 1);
				if (e < be) {
					be = e;
					best_decomp = {TState(br - 1, up - 2, i, k), TState(k + 2, j - 1)};
				}
			}
		}

		// Coaxial stack decomposition.
		if (br >= 2) {
			e = ML[br - 2][up][i][k] + CxFl[k + 1][j];
			if (e < be) {
				be = e;
				best_decomp = {TState(br - 2, up, i, k), TState(CxFlT, k + 1, j)};
			}
			if (up >= 2) {
				e = ML[br - 2][up - 2][i][k] + CxMM[k + 1][j];
				if (e < be) {
					be = e;
					best_decomp = {TState(br - 2, up - 2, i, k), TState(CxMMT, k + 1, j)};
				}
			}
		}

	}

	for (const TState &state : best_decomp)
		s.push(state);
}

librnary::Matching librnary::UpBrFolder::Traceback() {
	const int N = static_cast<int>(rna.size());
	auto m = EmptyMatching(static_cast<unsigned>(rna.size()));
	if (N == 0)
		return m;
	for (int i = 0; i < N; ++i)
		m[i] = i;
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
		else if (s.top().t == CxMMT)
			TraceCxMM(s);
		else // s.top().t == MLT
			TraceML(s);
	}
	return m;
}


int librnary::UpBrFolder::Fold(const PrimeStructure &primary) {
	// Init sequence data.
	rna = primary;
	em.SetRNA(primary);

	// The length of the RNA.
	const unsigned long RSZ = rna.size();
	// N is a convenience variable.
	// I have it as well as RSZ to avoid type warnings.
	// At least forcing this cast will catch any overflow errors, potentially.
	const int N = static_cast<int>(RSZ);

	if (N == 0)
		return 0;

	// Clear the DP tables.
	P.clear();
	ML.clear();
	E.clear();
	CxFl.clear();
	CxMM.clear();

	// Reset the DP tables.
	P.resize(RSZ, VE(RSZ, em.MaxMFE()));
	ML.resize((unsigned long) min(N / 2, max_multi_branches) + 1,
			  VVVE((unsigned long) min(N, max_multi_unpaired) + 1,
				   VVE(RSZ, VE(RSZ, em.MaxMFE()))));
	E.resize(RSZ, 0);
	CxFl.resize(RSZ, VE(RSZ, em.MaxMFE()));
	CxMM.resize(RSZ, VE(RSZ, em.MaxMFE()));

	// Special base case for ML table.
	// Can only end on a single unpaired nucleotide. No other states are base cases.
	for (int i = 0; i < N; ++i)
		ML[0][1][i][i] = 0;

	for (int i = N - 2; i >= 0; --i) { // i is 5' nucleotide.
		for (int j = i + 1; j < N; ++j) { // j is 3' nucleotide.
			if (ValidPair(rna[i], rna[j]) &&
				(lonely_pairs || !MustBeLonelyPair(rna, i, j, em.MIN_HAIRPIN_UNPAIRED))) {
				energy_t best = em.OneLoop(i, j); // Hairpin.
				for (int up = 0; up < min(max_multi_unpaired + 1, j - i); ++up) { // Multiloops.
					for (int br = 2; br <= min(max_multi_branches, (j - i - 1 - up) / 2); ++br) {
						int init = em.MLClosure(up, br + 1) + em.Branch(i, j);
						best = min(best, ML[br][up][i + 1][j - 1] + init);
						if (i + 2 < j - 1 && up >= 1) // Left dangle.
							best = min(best, ML[br][up - 1][i + 2][j - 1]
								+ init + em.ClosingThreeDangle(i, j));
						if (i + 1 < j - 2 && up >= 1) // Right dangle.
							best = min(best, ML[br][up - 1][i + 1][j - 2]
								+ init + em.ClosingFiveDangle(i, j));
						if (i + 2 < j - 2 && up >= 2) // Mismatch.
							best = min(best, ML[br][up - 2][i + 2][j - 2]
								+ init + em.ClosingMismatch(i, j));
						// Coaxial stack.
						for (int k = i + 1; k < j; ++k) {
							// Five prime.
							// ((_)_)
							//    ^ <- k
							if (k + 1 < j - 1 && i + 1 < k)
								best = min(best, ML[br - 1][up][k + 1][j - 1] + init
									+ em.FlushCoax(i, j, i + 1, k) +
									SSScore(i + 1, k));
							if (up >= 2) {
								// (.(_)_.)
								if (i + 2 < k && k + 1 < j - 2)
									best = min(best, ML[br - 1][up - 2][k + 1][j - 2] + init
										+
											em.MismatchCoax(i, j, i + 2, k) +
										SSScore(i + 2, k));
								// (.(_)._)
								if (i + 2 < k && k + 2 < j)
									best = min(best, ML[br - 1][up - 2][k + 2][j - 1] + init
										+
											em.MismatchCoax(i + 2, k, i, j) +
										SSScore(i + 2, k));
							}
							// Three prime.
							// (_(_))
							//   ^ <- k
							if (i + 1 < k - 1 && k < j - 1)
								best = min(best, ML[br - 1][up][i + 1][k - 1] + init
									+ em.FlushCoax(i, j, k, j - 1) +
									SSScore(k, j - 1));
							if (up >= 2) {
								// (._(_).)
								if (k < j - 2 && i + 2 < k - 1)
									best = min(best, ML[br - 1][up - 2][i + 2][k - 1] + init
										+
											em.MismatchCoax(i, j, k, j - 2) +
										SSScore(k, j - 2));
								// (_.(_).)
								if (k < j - 2 && i + 1 < k - 2)
									best = min(best, ML[br - 1][up - 2][i + 1][k - 2] + init
										+
											em.MismatchCoax(k, j - 2, i, j) +
										SSScore(k, j - 2));
							}
						}
					}
				}
				// Two loops for the two-loops (stack, bulge, or internal loop).
				for (int k = i + 1; k + 1 < j && (k - i - 1) <= max_twoloop_unpaired; ++k)
					for (int l = j - 1; l > k && (j - l - 1) + (k - i - 1) <= max_twoloop_unpaired; --l)
						best = min(best, P[k][l] + em.TwoLoop(i, k, l, j));
				P[i][j] = best;
			}

			// Fill the coaxial stack tables.
			for (int k = i + 1; k + 1 < j; ++k) {
				CxFl[i][j] = min(CxFl[i][j], em.FlushCoax(i, k, k + 1, j)
					+ SSScore(i, k) + SSScore(k + 1, j));
				if (i + 1 < k - 1)
					CxMM[i][j] = min(CxMM[i][j], em.MismatchCoax(i + 1, k - 1, k + 1, j)
						+ SSScore(i + 1, k - 1) + SSScore(k + 1, j));
				if (k + 2 < j - 1)
					CxMM[i][j] = min(CxMM[i][j], em.MismatchCoax(k + 2, j - 1, i, k)
						+ SSScore(i, k) + SSScore(k + 2, j - 1));
			}

			for (int up = 0;
				 up <= min(max_multi_unpaired, j - i + 1); ++up) { // up is the number of unpaired needed for valid ML.
				for (int br = 0; br <= min(max_multi_branches, (j - i + 1 - up) / 2); ++br) {
					energy_t best = em.MaxMFE();
					if (up > 0) // Single stranded.
						best = ML[br][up - 1][i][j - 1];
					if (br == 1) { // End on branch cases.
						if (up == 0)
							best = min(best, SSScore(i, j));
						if (i + 1 < j && up == 1) {
							best = min(best,
									   SSScore(i + 1, j) + em.FiveDangle(i + 1, j));
							best = min(best,
									   SSScore(i, j - 1) + em.ThreeDangle(i, j - 1));
						}
						if (i + 1 < j - 1 && up == 2) {
							best = min(best, SSScore(i + 1, j - 1) +
								em.Mismatch(i + 1, j - 1));
						}
					}
					// End on coaxial stack.
					if (br == 2) {
						if (up == 2)
							best = min(best, CxMM[i][j]);
						else if (up == 0)
							best = min(best, CxFl[i][j]);
					}
					for (int k = i; k + 2 <= j; ++k) { // Try all decompositions into 5' ML fragment and 3' branch.
						if (br >= 1) {
							best = min(best, ML[br - 1][up][i][k] + SSScore(k + 1, j));
							if (k + 2 < j && up >= 1)
								best = min(best, ML[br - 1][up - 1][i][k] + SSScore(k + 2, j) +
									em.FiveDangle(k + 2, j));
							if (k + 1 < j - 1 && up >= 1)
								best = min(best, ML[br - 1][up - 1][i][k] + SSScore(k + 1, j - 1) +
									em.ThreeDangle(k + 1, j - 1));
							if (k + 2 < j - 1 && up >= 2)
								best = min(best, ML[br - 1][up - 2][i][k] + SSScore(k + 2, j - 1) +
									em.Mismatch(k + 2, j - 1));
						}
						if (br >= 2) {
							// Coaxial stack decomposition.
							best = min(best, ML[br - 2][up][i][k] + CxFl[k + 1][j]);
							if (up >= 2)
								best = min(best, ML[br - 2][up - 2][i][k] + CxMM[k + 1][j]);
						}
					}
					ML[br][up][i][j] = best;
				}
			}
		}
	}

	for (int i = 1; i < N; ++i) {
		energy_t best = E[i - 1];
		for (int k = -1; k < i; ++k) {
			energy_t decomp = k == -1 ? 0 : E[k];
			best = min(best, decomp + SSScore(k + 1, i));
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
		E[i] = best;
	}

	return E[N - 1];
}

librnary::VVE librnary::UpBrFolder::GetP() const {
	return P;
}

librnary::VE librnary::UpBrFolder::GetE() const {
	return E;
}

librnary::StrainedUnpairedModel librnary::UpBrFolder::GetEM() const {
	return em;
}

void librnary::UpBrFolder::SetMaxBranches(unsigned v) {
	max_multi_branches = v;
}

int librnary::UpBrFolder::MaxBranches() const {
	return max_multi_branches;
}

void librnary::UpBrFolder::SetMaxMulti(unsigned max_up) {
	max_multi_unpaired = max_up;
}

int librnary::UpBrFolder::MaxMulti() const {
	return max_multi_unpaired;
}

void librnary::UpBrFolder::SetMaxTwoLoop(unsigned max_up) {
	max_twoloop_unpaired = max_up;
}

int librnary::UpBrFolder::MaxTwoLoop() const {
	return max_twoloop_unpaired;
}

void librnary::UpBrFolder::SetLonelyPairs(bool v) {
	lonely_pairs = v;
}

bool librnary::UpBrFolder::LonelyPairs() const {
	return lonely_pairs;
}

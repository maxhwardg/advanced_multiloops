//
// Created by max on 6/28/16.
//

#include "folders/aalberts_folder.hpp"

using namespace std;

librnary::energy_t librnary::AalbertsFolder::SSScore(int i, int j) {
	return em.Branch(i, j) + P[i][j];
}

void librnary::AalbertsFolder::TraceE(stack<TState> &s) {
	// Manage current state.
	TState curr = s.top();
	int i = curr.i;
	s.pop();

	// Base case.
	if (i == 0)
		return;

	// Maintain best decomposition.
	vector<TState> best_decomp = {TState(i - 1)}; // Unpaired 3'.
	librnary::energy_t be = E[i - 1];
	for (int k = -1; k < i; ++k) {
		int decomp = k == -1 ? 0 : E[k];
		vector<TState> decomp_state;
		if (k != -1)
			decomp_state.push_back(TState(k));
		librnary::energy_t e = decomp + SSScore(k + 1, i);
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
		e = decomp + Cx[k + 1][i];
		if (e < be) {
			be = e;
			best_decomp = {TState(CxT, k + 1, i)};
			best_decomp.insert(best_decomp.end(), decomp_state.begin(), decomp_state.end());
		}

	}

	for (const TState &state : best_decomp)
		s.push(state);
}

void librnary::AalbertsFolder::TraceP(stack<TState> &s) {
	TState curr = s.top();
	s.pop();
	int i = curr.i, j = curr.j;

	// Hairpin. Also technically a base case.
	vector<TState> best_decomp; // Unpaired 3'.
	librnary::energy_t be = em.OneLoop(i, j);
	librnary::energy_t e;
	// Multi-loops.
	for (int a = 0; a <= MaxALengthSegs(j - i - 1); ++a) {
		for (int b = 0; b <= MaxBLengthSegs(j - i - 1); ++b) {
			int init_nocx = em.MLInit(a + 1, b + 1) + em.Branch(i, j);
			int init_cx = em.MLInit(a + 2, b) + em.Branch(i, j);
			// Vanilla.
			e = ML[2][b][a][i + 1][j - 1] + init_nocx;
			if (e < be) {
				be = e;
				best_decomp = {TState(2, b, a, i + 1, j - 1)};
			}
			if (i + 2 < j - 1 && a >= 1) { // Left dangle.
				e = ML[2][b][a - 1][i + 2][j - 1] + init_nocx + em.ClosingThreeDangle(i, j);
				if (e < be) {
					be = e;
					best_decomp = {TState(2, b, a - 1, i + 2, j - 1)};
				}
			}
			if (i + 1 < j - 2 && a >= 1) { // Right dangle.
				e = ML[2][b][a - 1][i + 1][j - 2] + init_nocx + em.ClosingFiveDangle(i, j);
				if (e < be) {
					be = e;
					best_decomp = {TState(2, b, a - 1, i + 1, j - 2)};
				}
			}
			if (i + 2 < j - 2 && a >= 2) { // Mismatch.
				e = ML[2][b][a - 2][i + 2][j - 2] + init_nocx + em.ClosingMismatch(i, j);
				if (e < be) {
					be = e;
					best_decomp = {TState(2, b, a - 2, i + 2, j - 2)};
				}
			}


			// Coaxial stack.
			for (int k = i + 1; k < j; ++k) {
				// Five prime.
				// ((_)_)
				//    ^ <- k
				if (k + 1 < j - 1 && i + 1 < k) {
					e = ML[1][b][a][k + 1][j - 1] + init_cx + em.FlushCoax(i, j, i + 1, k) + SSScore(i + 1, k);
					if (e < be) {
						be = e;
						best_decomp = {TState(1, b, a, k + 1, j - 1), TState(i + 1, k)};
					}
				}
				// (.(_)_.)
				if (i + 2 < k && k + 1 < j - 2) {
					e = ML[1][b][a][k + 1][j - 2] + init_cx + em.MismatchCoax(i, j, i + 2, k) + SSScore(i + 2, k);
					if (e < be) {
						be = e;
						best_decomp = {TState(1, b, a, k + 1, j - 2), TState(i + 2, k)};
					}
				}
				// (.(_)._)
				if (i + 2 < k && k + 2 < j) {
					e = ML[1][b][a][k + 2][j - 1] + init_cx + em.MismatchCoax(i + 2, k, i, j) + SSScore(i + 2, k);
					if (e < be) {
						be = e;
						best_decomp = {TState(1, b, a, k + 2, j - 1), TState(i + 2, k)};
					}
				}

				// Three prime.
				// (_(_))
				//   ^ <- k
				if (i + 1 < k - 1 && k < j - 1) {
					e = ML[1][b][a][i + 1][k - 1] + init_cx + em.FlushCoax(i, j, k, j - 1) + SSScore(k, j - 1);
					if (e < be) {
						be = e;
						best_decomp = {TState(1, b, a, i + 1, k - 1), TState(k, j - 1)};
					}
				}
				// (._(_).)
				if (k < j - 2 && i + 2 < k - 1) {
					e = ML[1][b][a][i + 2][k - 1] + init_cx + em.MismatchCoax(i, j, k, j - 2) + SSScore(k, j - 2);
					if (e < be) {
						be = e;
						best_decomp = {TState(1, b, a, i + 2, k - 1), TState(k, j - 2)};
					}
				}
				// (_.(_).)
				if (k < j - 2 && i + 1 < k - 2) {
					e = ML[1][b][a][i + 1][k - 2] + init_cx + em.MismatchCoax(k, j - 2, i, j) + SSScore(k, j - 2);
					if (e < be) {
						be = e;
						best_decomp = {TState(1, b, a, i + 1, k - 2), TState(k, j - 2)};
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

void librnary::AalbertsFolder::TraceCx(stack<TState> &s) {
	TState curr = s.top();
	s.pop();
	int i = curr.i, j = curr.j;

	vector<TState> best_decomp;
	librnary::energy_t be = em.MaxMFE();
	librnary::energy_t e;

	for (int k = i + 1; k + 1 < j; ++k) {
		e = em.FlushCoax(i, k, k + 1, j) + SSScore(i, k) + SSScore(k + 1, j);
		if (e < be) {
			be = e;
			best_decomp = {TState(i, k), TState(k + 1, j)};
		}
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


void librnary::AalbertsFolder::TraceML(stack<TState> &s) {
	TState curr = s.top();
	s.pop();
	int br = curr.br, b = curr.b, a = curr.a, i = curr.i, j = curr.j;

	// Base case.
	if (i == j)
		return;

	vector<TState> best_decomp;
	librnary::energy_t be = em.MaxMFE();

	if (a >= 1) { // Single stranded.
		best_decomp = {TState(br, b, a - 1, i, j - 1)};
		be = ML[br][b][a - 1][i][j - 1];
	}

	librnary::energy_t e;
	if (br <= 1 && b == 1) { // End on branch cases.
		if (a == 1) {
			e = SSScore(i, j);
			if (e < be) {
				be = e;
				best_decomp = {TState(i, j)};
			}
		}
		if (i + 1 < j && a == 2) {
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
		if (i + 1 < j - 1 && a == 3) {
			e = SSScore(i + 1, j - 1) + em.Mismatch(i + 1, j - 1);
			if (e < be) {
				be = e;
				best_decomp = {TState(i + 1, j - 1)};
			}
		}
	}

	// End on coaxial stack.
	if (br <= 2 && a == 2 && b == 0) {
		e = Cx[i][j];
		if (e < be) {
			be = e;
			best_decomp = {TState(CxT, i, j)};
		}
	}

	for (int k = i; k + 2 <= j; ++k) { // Try all decompositions into 5' ML fragment and 3' branch.
		int brprime = max(0, br - 1);
		if (b >= 1) {
			if (a >= 1) {
				e = ML[brprime][b - 1][a - 1][i][k] + SSScore(k + 1, j);
				if (e < be) {
					be = e;
					best_decomp = {TState(brprime, b - 1, a - 1, i, k), TState(k + 1, j)};
				}
			}
			if (k + 2 < j && a >= 2) {
				e = ML[brprime][b - 1][a - 2][i][k] + SSScore(k + 2, j) + em.FiveDangle(k + 2, j);
				if (e < be) {
					be = e;
					best_decomp = {TState(brprime, b - 1, a - 2, i, k), TState(k + 2, j)};
				}
			}
			if (k + 1 < j - 1 && a >= 2) {
				e = ML[brprime][b - 1][a - 2][i][k] + SSScore(k + 1, j - 1) +
					em.ThreeDangle(k + 1, j - 1);
				if (e < be) {
					be = e;
					best_decomp = {TState(brprime, b - 1, a - 2, i, k), TState(k + 1, j - 1)};
				}
			}
			if (k + 2 < j - 1 && a >= 3) {
				e = ML[brprime][b - 1][a - 3][i][k] + SSScore(k + 2, j - 1) +
					em.Mismatch(k + 2, j - 1);
				if (e < be) {
					be = e;
					best_decomp = {TState(brprime, b - 1, a - 3, i, k), TState(k + 2, j - 1)};
				}
			}
		}

		// Coaxial stack decomposition.
		if (a >= 2) {
			e = ML[0][b][a - 2][i][k] + Cx[k + 1][j];
			if (e < be) {
				be = e;
				best_decomp = {TState(0, b, a - 2, i, k), TState(CxT, k + 1, j)};
			}
		}

	}

	for (const TState &state : best_decomp)
		s.push(state);
}

librnary::Matching librnary::AalbertsFolder::Traceback() {
	const int N = static_cast<int>(rna.size());
	auto m = EmptyMatching(static_cast<unsigned>(N));
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
		} else if (s.top().t == CxT)
			TraceCx(s);
		else // s.top().t == MLT
			TraceML(s);
	}
	return m;
}

int librnary::AalbertsFolder::Fold(const PrimeStructure &primary) {
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
	Cx.clear();

	// Reset the DP tables.
	P.resize(RSZ, V<librnary::energy_t>(RSZ, em.MaxMFE()));
	ML.resize(3, _4DVE((unsigned long) MaxBLengthSegs(N) + 1,
					   VVVE((unsigned long) MaxALengthSegs(N) + 1, VVE(RSZ, VE(RSZ, em.MaxMFE())))));
	E.resize(RSZ, 0);
	Cx.resize(RSZ, VE(RSZ, em.MaxMFE()));

	// Special base case for ML table.
	// Can only end on a single unpaired nucleotide. No other states are base cases.
	for (int i = 0; i < N; ++i)
		ML[0][0][1][i][i] = 0;

	for (int i = N - 2; i >= 0; --i) { // i is 5' nucleotide.
		for (int j = i + 1; j < N; ++j) { // j is 3' nucleotide.
			if (ValidPair(rna[i], rna[j]) &&
				(lonely_pairs || !MustBeLonelyPair(rna, i, j, em.MIN_HAIRPIN_UNPAIRED))) {
				librnary::energy_t best = em.OneLoop(i, j); // Hairpins.
				// Multiloops.
				for (int a = 0; a <= MaxALengthSegs(j - i - 1); ++a) {
					for (int b = 0; b <= MaxBLengthSegs(j - i - 1); ++b) {
						// The +1s and +2s here account for the closing branch.
						int init_nocx = em.MLInit(a + 1, b + 1) + em.Branch(i, j);
						int init_cx = em.MLInit(a + 2, b) + em.Branch(i, j);
						// (_)
						best = min(best, ML[2][b][a][i + 1][j - 1] + init_nocx);
						if (i + 2 < j - 1 && a >= 1) // Left dangle.
							best = min(best, ML[2][b][a - 1][i + 2][j - 1]
								+ init_nocx + em.ClosingThreeDangle(i, j));
						if (i + 1 < j - 2 && a >= 1) // Right dangle.
							best = min(best, ML[2][b][a - 1][i + 1][j - 2]
								+ init_nocx + em.ClosingFiveDangle(i, j));
						if (i + 2 < j - 2 && a >= 2) // Mismatch.
							best = min(best, ML[2][b][a - 2][i + 2][j - 2]
								+ init_nocx + em.ClosingMismatch(i, j));
						// Coaxial stack.
						for (int k = i + 1; k < j; ++k) {
							// Five prime.
							// ((_)_)
							//    ^ <- k
							if (k + 1 < j - 1 && i + 1 < k)
								best = min(best,
										   ML[1][b][a][k + 1][j - 1] + init_cx + em.FlushCoax(i, j, i + 1, k)
											   + SSScore(i + 1, k));
							// (.(_)_.)
							if (i + 2 < k && k + 1 < j - 2)
								best = min(best,
										   ML[1][b][a][k + 1][j - 2] + init_cx + em.MismatchCoax(i, j, i + 2, k)
											   + SSScore(i + 2, k));
							// (.(_)._)
							if (i + 2 < k && k + 2 < j)
								best = min(best,
										   ML[1][b][a][k + 2][j - 1] + init_cx + em.MismatchCoax(i + 2, k, i, j)
											   + SSScore(i + 2, k));
							// Three prime.
							// (_(_))
							//   ^ <- k
							if (i + 1 < k - 1 && k < j - 1)
								best = min(best,
										   ML[1][b][a][i + 1][k - 1] + init_cx + em.FlushCoax(i, j, k, j - 1)
											   + SSScore(k, j - 1));
							// (._(_).)
							if (k < j - 2 && i + 2 < k - 1)
								best = min(best,
										   ML[1][b][a][i + 2][k - 1] + init_cx + em.MismatchCoax(i, j, k, j - 2)
											   + SSScore(k, j - 2));
							// (_.(_).)
							if (k < j - 2 && i + 1 < k - 2)
								best = min(best,
										   ML[1][b][a][i + 1][k - 2] + init_cx + em.MismatchCoax(k, j - 2, i, j)
											   + SSScore(k, j - 2));
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
				Cx[i][j] = min(Cx[i][j], em.FlushCoax(i, k, k + 1, j)
					+ SSScore(i, k) + SSScore(k + 1, j));
				if (i + 1 < k - 1)
					Cx[i][j] = min(Cx[i][j], em.MismatchCoax(i + 1, k - 1, k + 1, j)
						+ SSScore(i + 1, k - 1) + SSScore(k + 1, j));
				if (k + 2 < j - 1)
					Cx[i][j] = min(Cx[i][j], em.MismatchCoax(k + 2, j - 1, i, k)
						+ SSScore(i, k) + SSScore(k + 2, j - 1));
			}

			for (int br = 0; br <= 2; ++br) {
				for (int a = 0; a <= MaxALengthSegs(j - i + 1); ++a) {
					for (int b = 0; b <= MaxBLengthSegs(j - i + 1); ++b) {
						int best = em.MaxMFE();
						if (a >= 1) // Single stranded.
							best = ML[br][b][a - 1][i][j - 1];
						if (br <= 1 && b == 1) { // End on branch cases.
							// Vanilla.
							if (a == 1)
								best = min(best, SSScore(i, j));
							// Dangles.
							if (i + 1 < j && a == 2) {
								best = min(best, SSScore(i + 1, j) + em.FiveDangle(i + 1, j));
								best = min(best, SSScore(i, j - 1) + em.ThreeDangle(i, j - 1));
							}
							// Mismatch
							if (i + 1 < j - 1 && a == 3) {
								best = min(best, SSScore(i + 1, j - 1) + em.Mismatch(i + 1, j - 1));
							}
						}
						// End on coaxial stack.
						if (br <= 2 && a == 2 && b == 0) {
							best = min(best, Cx[i][j]);
						}
						// Try all decompositions into 5' ML fragment and 3' branch.
						for (int k = i; k + 2 <= j; ++k) {
							int brprime = max(0, br - 1);
							// Non-coax decomps.
							if (b >= 1) {
								// {_}(_)
								if (a >= 1)
									best = min(best, ML[brprime][b - 1][a - 1][i][k] + SSScore(k + 1, j));
								// {_}.(_)
								if (k + 2 < j && a >= 2)
									best = min(best,
											   ML[brprime][b - 1][a - 2][i][k] + SSScore(k + 2, j)
												   + em.FiveDangle(k + 2, j));
								// {_}(_).
								if (k + 1 < j - 1 && a >= 2)
									best = min(best,
											   ML[brprime][b - 1][a - 2][i][k] + SSScore(k + 1, j - 1)
												   + em.ThreeDangle(k + 1, j - 1));
								// {_}.(_).
								if (k + 2 < j - 1 && a >= 3)
									best = min(best,
											   ML[brprime][b - 1][a - 3][i][k] + SSScore(k + 2, j - 1)
												   + em.Mismatch(k + 2, j - 1));
							}
							if (a >= 2) {
								// Coaxial stack decomposition.
								best = min(best, ML[0][b][a - 2][i][k] + Cx[k + 1][j]);
							}
						}
						ML[br][b][a][i][j] = best;
					}
				}
			}
		}
	}

	for (int i = 1; i < N; ++i) {
		librnary::energy_t best = E[i - 1];
		for (int k = -1; k < i; ++k) {
			int decomp = k == -1 ? 0 : E[k];
			best = min(best, decomp + SSScore(k + 1, i));
			if (k + 2 < i)
				best = min(best, decomp + SSScore(k + 2, i) + em.FiveDangle(k + 2, i));
			if (k + 1 < i - 1)
				best = min(best, decomp + SSScore(k + 1, i - 1) + em.ThreeDangle(k + 1, i - 1));
			if (k + 2 < i - 1)
				best = min(best, decomp + SSScore(k + 2, i - 1) + em.Mismatch(k + 2, i - 1));
			// Coaxial stack decompositions.
			best = min(best, decomp + Cx[k + 1][i]);
		}
		E[i] = best;
	}

	return E[N - 1];
}

librnary::VVE librnary::AalbertsFolder::GetP() const {
	return P;
}

librnary::VE librnary::AalbertsFolder::GetE() const {
	return E;
}

librnary::AalbertsModel librnary::AalbertsFolder::GetEM() const {
	return em;
}

void librnary::AalbertsFolder::SetMaxBLength(unsigned v) {
	max_multi_blength = v;
}

int librnary::AalbertsFolder::MaxBLength() const {
	return max_multi_blength;
}

void librnary::AalbertsFolder::SetMaxALength(unsigned max_up) {
	max_multi_alength = max_up;
}

int librnary::AalbertsFolder::MaxALength() const {
	return max_multi_alength;
}

void librnary::AalbertsFolder::SetMaxTwoLoop(unsigned max_up) {
	max_twoloop_unpaired = max_up;
}

int librnary::AalbertsFolder::MaxTwoLoop() const {
	return max_twoloop_unpaired;
}

void librnary::AalbertsFolder::SetLonelyPairs(bool v) {
	lonely_pairs = v;
}

bool librnary::AalbertsFolder::LonelyPairs() const {
	return lonely_pairs;
}
void librnary::AalbertsFolder::SetModel(const librnary::AalbertsModel &_em) {
	this->em = _em;
}

//
// Created by max on 5/30/16.
//

#include "statistics.hpp"

#include <set>

using namespace std;

int librnary::slippage::TruePositives(const librnary::Matching &matching_true,
									  const librnary::Matching &matching_proband) {
	assert(matching_proband.size() == matching_true.size());

	const int N = static_cast<int>(matching_true.size());

	set<pair<int, int>> true_pairs;
	for (int i = 0; i < N; ++i) {
		if (i < matching_true[i]) {
			true_pairs.emplace(i, matching_true[i]);
		}
	}

	int tps = 0;

	for (int i = 0; i < N; ++i) {
		if (i < matching_proband[i]) {
			if (true_pairs.count({i, matching_proband[i]})
				|| true_pairs.count({i + 1, matching_proband[i]})
				|| true_pairs.count({i - 1, matching_proband[i]})
				|| true_pairs.count({i, matching_proband[i] + 1})
				|| true_pairs.count({i, matching_proband[i] - 1})) {
				++tps;
			}
		}
	}
	return tps;
}

int librnary::slippage::FalseNegatives(const librnary::Matching &matching_true,
									   const librnary::Matching &matching_proband) {
	assert(matching_proband.size() == matching_true.size());

	const int N = static_cast<int>(matching_true.size());


	set<pair<int, int>> proband_pairs;
	for (int i = 0; i < N; ++i) {
		if (i < matching_proband[i]) {
			proband_pairs.emplace(i, matching_proband[i]);
		}
	}

	int fns = 0;

	for (int i = 0; i < N; ++i) {
		if (i < matching_true[i]) {
			if (!proband_pairs.count({i, matching_true[i]})
				&& !proband_pairs.count({i + 1, matching_true[i]})
				&& !proband_pairs.count({i - 1, matching_true[i]})
				&& !proband_pairs.count({i, matching_true[i] + 1})
				&& !proband_pairs.count({i, matching_true[i] - 1})) {
				++fns;
			}
		}
	}
	return fns;

}

int librnary::slippage::FalsePositives(const librnary::Matching &matching_true,
									   const librnary::Matching &matching_proband) {
	return FalseNegatives(matching_proband, matching_true);
}

double librnary::slippage::Sensitivity(const librnary::Matching &mtrue, const librnary::Matching &mpred) {
	double tp = TruePositives(mtrue, mpred);
	double fn = FalseNegatives(mtrue, mpred);
	if (tp <= 0.0) return 0.0;
	return tp / (tp + fn);
}

double librnary::slippage::PositivePredictiveValue(const librnary::Matching &mtrue, const librnary::Matching &mpred) {
	double tp = TruePositives(mtrue, mpred);
	double fp = FalsePositives(mtrue, mpred);
	if (tp <= 0.0) return 0.0;
	return tp / (tp + fp);
}

double librnary::slippage::F1Score(const librnary::Matching &mtrue, const librnary::Matching &mpred) {
	double tp = TruePositives(mtrue, mpred);
	double fp = FalsePositives(mtrue, mpred);
	double fn = FalseNegatives(mtrue, mpred);
	return (2 * tp) / (2 * tp + fp + fn);
}

int librnary::TruePositives(const librnary::Matching &matching_true, const librnary::Matching &matching_proband) {
	assert(matching_proband.size() == matching_true.size());
	int tps = 0;
	std::vector<bool> marked(matching_true.size(), false);

	for (int i = 0; i < static_cast<int>(matching_true.size()); ++i) {
		if (marked[i])
			continue;
		if (matching_true[i] != i && matching_true[i] == matching_proband[i]) {
			tps += 1;
			marked[matching_true[i]] = true;
		}
	}
	return tps;
}

int librnary::FalseNegatives(const librnary::Matching &matching_true, const librnary::Matching &matching_proband) {
	assert(matching_proband.size() == matching_true.size());
	int fns = 0;
	std::vector<bool> marked(matching_true.size(), false);

	for (int i = 0; i < static_cast<int>(matching_true.size()); ++i) {
		if (marked[i])
			continue;
		if (matching_true[i] != i && matching_proband[i] != matching_true[i]) {
			fns += 1;
			marked[matching_true[i]] = true;
		}
	}
	return fns;
}

int librnary::FalsePositives(const librnary::Matching &matching_true, const librnary::Matching &matching_proband) {
	return FalseNegatives(matching_proband, matching_true);
}

double librnary::Sensitivity(const librnary::Matching &mtrue, const librnary::Matching &mpred) {
	double tp = TruePositives(mtrue, mpred);
	double fn = FalseNegatives(mtrue, mpred);
	if (tp <= 0.0) return 0.0;
	return tp / (tp + fn);
}

double librnary::PositivePredictiveValue(const librnary::Matching &mtrue, const librnary::Matching &mpred) {
	double tp = TruePositives(mtrue, mpred);
	double fp = FalsePositives(mtrue, mpred);
	if (tp <= 0.0) return 0.0;
	return tp / (tp + fp);
}

double librnary::F1Score(const librnary::Matching &mtrue, const librnary::Matching &mpred) {
	double tp = TruePositives(mtrue, mpred);
	double fp = FalsePositives(mtrue, mpred);
	double fn = FalseNegatives(mtrue, mpred);
	return (2 * tp) / (2 * tp + fp + fn);
}

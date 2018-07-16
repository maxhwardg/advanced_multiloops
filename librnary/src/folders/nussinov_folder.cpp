//
// Created by max on 8/8/16.
//

#include "folders/nussinov_folder.hpp"

#include <stack>

using namespace std;

int librnary::NussinovFolder::Fold(size_t N, std::function<int(int, int)> _score_fn) {
	this->score_fn = _score_fn;
	RNA_LEN = static_cast<int>(N);
	if (RNA_LEN == 0)
		return 0;
	F.assign(static_cast<size_t>(RNA_LEN + 1), VI(static_cast<size_t>(RNA_LEN), 0));
	for (int i = RNA_LEN - 1; i >= 0; --i) {
		for (int j = i + 1; j < RNA_LEN; ++j) {
			int best = score_fn(i, j) + F[i + 1][j - 1];
			for (int k = i; k < j; ++k) {
				best = max(best, F[i][k] + F[k + 1][j]);
			}
			F[i][j] = best;
		}
	}
	return F[0][RNA_LEN - 1];
}

librnary::Matching librnary::NussinovFolder::Traceback() {
	if (RNA_LEN == 0)
		return EmptyMatching(0);
	stack<pair<int, int>> s;
	s.push({0, RNA_LEN - 1});
	librnary::Matching match = EmptyMatching(static_cast<unsigned>(RNA_LEN));
	while (!s.empty()) {
		int i = s.top().first, j = s.top().second;
		s.pop();
		if (i >= j)
			continue;
		int best = score_fn(i, j) + F[i + 1][j - 1];
		vector<pair<int, int>> decomp = {{i + 1, j - 1}};
		for (int k = i; k < j; ++k) {
			int tmp_sc = F[i][k] + F[k + 1][j];
			if (tmp_sc > best) {
				best = tmp_sc;
				decomp = {{i, k}, {k + 1, j}};
			}
		}
		for (const auto &p : decomp) {
			s.push(p);
		}
		if (decomp.size() == 1) {
			match[i] = j;
			match[j] = i;
		}
	}
	return match;
}
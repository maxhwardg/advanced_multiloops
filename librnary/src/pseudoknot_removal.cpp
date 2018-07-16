//
// Created by max on 8/9/16.
//

#include <folders/nussinov_folder.hpp>
#include "pseudoknot_removal.hpp"

using namespace std;

librnary::Matching librnary::RemovePseudoknotsMaximizePairs(const librnary::Matching &match) {
	librnary::NussinovFolder folder;
	folder.Fold(match.size(), [&](int i, int j) -> int {
		return match[i] == j ? 1 : -1;
	});
	return folder.Traceback();
}
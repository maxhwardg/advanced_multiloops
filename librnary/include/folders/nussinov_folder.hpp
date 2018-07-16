//
// Created by max on 8/8/16.
//

#ifndef RNARK_NUSSINOV_FOLDER_HPP
#define RNARK_NUSSINOV_FOLDER_HPP

#include "primary_structure.hpp"
#include "secondary_structure.hpp"
#include "vector_types.hpp"

#include <vector>
#include <functional>

namespace librnary {
/**
 * An implementation of a Nussinov base pair maximization algorithm.
 */

class NussinovFolder {
	int RNA_LEN;
	std::function<int(int, int)> score_fn;

	/// Best substructure in the fragment (hence F) [i,j].
	VVI F;
public:
	int Fold(size_t N, std::function<int(int, int)> score_fn);
	Matching Traceback();
};

}

#endif //RNARK_NUSSINOV_FOLDER_HPP

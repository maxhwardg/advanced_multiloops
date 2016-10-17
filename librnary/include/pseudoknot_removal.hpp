//
// Created by max on 8/9/16.
//

#ifndef RNARK_PSEUDOKNOT_REMOVAL_HPP
#define RNARK_PSEUDOKNOT_REMOVAL_HPP

#include "secondary_structure.hpp"
namespace librnary {

/**
 * Removes pseudoknots from a matching in a way that results in the maximum remaining pairs.
 * @param match The matching.
 * @return A pseudoknot free matching with maixmum reamining pairs.
 */
Matching RemovePseudoknotsMaximizePairs(const Matching &match);

}

#endif //RNARK_PSEUDOKNOT_REMOVAL_HPP

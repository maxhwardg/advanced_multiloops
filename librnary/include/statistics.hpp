//
// Created by max on 5/30/16.
//

/*
 * Contains functions for computing statistics relevant to RNA structure prediction.
 */

#ifndef RNARK_STATISTICS_HPP
#define RNARK_STATISTICS_HPP

#include "secondary_structure.hpp"

namespace librnary {

/// The slippage namespace has statistical functions that account for slippage.
/// See "Expanded Sequence Dependence of Thermodynamic Parameters Improves Prediction of RNA Secondary Structure"
/// For a defintion of slippage, and some reasons for considering it.
namespace slippage {
int TruePositives(const Matching &true_matching, const Matching &proband_matching);

int FalseNegatives(const Matching &true_matching, const Matching &proband_matching);

int FalsePositives(const Matching &true_matching, const Matching &proband_matching);

double Sensitivity(const Matching &true_matching, const Matching &proband_matching);

double PositivePredictiveValue(const Matching &true_matching, const Matching &proband_matching);

double F1Score(const Matching &true_matching, const Matching &proband_matching);
}

int TruePositives(const Matching &true_matching, const Matching &proband_matching);

int FalseNegatives(const Matching &true_matching, const Matching &proband_matching);

int FalsePositives(const Matching &true_matching, const Matching &proband_matching);

double Sensitivity(const Matching &true_matching, const Matching &proband_matching);

double PositivePredictiveValue(const Matching &true_matching, const Matching &proband_matching);

double F1Score(const Matching &true_matching, const Matching &proband_matching);
}

#endif //RNARK_STATISTICS_HPP

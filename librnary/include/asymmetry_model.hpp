//
// Created by max on 6/07/16.
//

#ifndef RNARK_ASYMMETRY_MODEL_HPP
#define RNARK_ASYMMETRY_MODEL_HPP

#include "nn_model.hpp"

namespace librnary {
/**
 * Computes multi-loop energy using asymmetry in unpaired nucleotides to either side of a branch as a parameter.
 * Implements an affine cost for asymmetry.
 * By default has no asymmetry code.
 */
class AsymmetryModel: public NNModel {
	energy_t ml_init = 93, ml_branch = -6, ml_unpaired = 0, ml_asymmetry = 0;
public:
	int SumAsymmetry(const librnary::Surface &surf) const;
	energy_t MLClosure(const librnary::SSTree &sst, SSTreeNodeId id) const override;
	energy_t MLClosure(const librnary::Surface &surf) const override;
	energy_t MLClosure(int branches, int unpaired, int asymmetry) const;
	energy_t MLInit() const;
	energy_t MLBranchCost() const;
	energy_t MLUnpairedCost() const;
	energy_t MLAsymmetryCost() const;
	void SetMLInit(energy_t v);
	void SetMLBranchCost(energy_t v);
	void SetMLUnpairedCost(energy_t v);
	void SetMLAsymmetryCost(energy_t v);
	void SetMLParams(energy_t ml_init,
					 energy_t ml_branch,
					 energy_t ml_unpaired,
					 energy_t ml_asymmetry);
	AsymmetryModel(const std::string &_data_path, const PrimeStructure &_rna)
		: NNModel(_data_path, _rna) {}
	AsymmetryModel(const std::string &_data_path)
		: NNModel(_data_path) {}
};
}

#endif //RNARK_ASYMMETRY_MODEL_HPP

//
// Created by max on 6/15/16.
//

#ifndef RNARK_NN_AFFINE_MODEL_HPP
#define RNARK_NN_AFFINE_MODEL_HPP

#include "nn_model.hpp"
#include "energy.hpp"

namespace librnary {
/**
 * Implements a typical nearest neighbour model with an affine function ofr multi-loop free energy.
 * Because RNAstructure is used as a back-end, the parameters from RNAstructure are used for multi-loops here.
 */
class NNAffineModel: public NNModel {
	/**
	 * Energy penalty for multi-loop initiation.
	 */
	energy_t ml_init;
	/**
	 * Energy penalty for multi-loop branch.
	 */
	energy_t ml_branch;
	/**
	 * Energy penalty for multi-loop unpaired.
	 */
	energy_t ml_unpiared;
public:
	energy_t MLClosure(const librnary::SSTree &sstree, librnary::SSTreeNodeId node_id) const override;
	energy_t MLClosure(const librnary::Surface &surf) const override;
	energy_t MLClosure(int branches, int unpaired) const;
	energy_t MLInitCost() const;
	virtual energy_t MLBranchCost() const;
	energy_t MLUnpairedCost() const;
	void SetMLInitCost(energy_t v);
	void SetMLBranchCost(energy_t v);
	void SetMLUnpairedCost(energy_t v);
	void SetMLParams(energy_t init, energy_t branch, energy_t unpaired);
	NNAffineModel(const std::string &data_path, const PrimeStructure &_rna)
		: NNModel(data_path, _rna) {
		ml_init = dt->efn2a;
		ml_branch = dt->efn2c;
		ml_unpiared = dt->efn2b;
	}
	NNAffineModel(const std::string &data_path)
		: NNModel(data_path) {
		ml_init = dt->efn2a;
		ml_branch = dt->efn2c;
		ml_unpiared = dt->efn2b;
	}
};
}

#endif //RNARK_NN_AFFINE_MODEL_HPP

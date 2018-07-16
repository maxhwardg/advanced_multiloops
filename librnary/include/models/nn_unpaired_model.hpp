//
// Created by max on 6/19/16.
//

#ifndef RNARK_NN_UNPAIRED_MODEL_HPP
#define RNARK_NN_UNPAIRED_MODEL_HPP

#include "nn_model.hpp"

namespace librnary {
/**
 * A nearest neighbour based energy model that allows arbitrarily complex energy functions involving free base pairs
 * in multi-branch loops. This, by default, is the classic logarithmic model as described in the Turner 1999 rules.
 */
class NNUnpairedModel: public NNModel {
protected:
	energy_t ml_init = 101, ml_br_cost = -3, ml_up_cost = -3;
	kcalmol_t ml_log_mult = 1.1;
	int ml_up_pivot = 6;
public:

	/// The constant cost of a branch in a multi-loop.
	energy_t MLBranchCost() const;
	/// The initiation cost of a multi-loop exlcluding the cost of branches.
	energy_t MLClosure(int unpaired) const;
	/// The initiation cost of a multi-loop given only the relevant features.
	virtual energy_t MLClosure(int unpaired, int branches) const;
	/// The multi-loop initiation cost of a given loop region. Assumed to be a multi-loop.
	virtual energy_t MLClosure(const librnary::SSTree &sstree, librnary::SSTreeNodeId node_id) const override;
	/// The multi-loop initiation cost of a given loop region. Assumed to be a multi-loop.
	virtual energy_t MLClosure(const librnary::Surface &surf) const override;

	void
	SetMLParams(energy_t ml_init, energy_t ml_br_cost, energy_t ml_up_cost, kcalmol_t ml_log_mult, int ml_up_pivot);
	energy_t MLInitConstant() const;
	void SetMLInitConstant(energy_t ml_init);
	void SetMLBranchCost(energy_t ml_br_cost);
	energy_t MLUnpairedCost() const;
	void SetMLUnpairedCost(energy_t ml_up_cost);
	kcalmol_t MLLogMultiplier() const;
	void SetMLLogMultiplier(kcalmol_t ml_log_mult);
	int MLUnpairedPivot() const;
	void SetMLUnpairedPivot(int ml_up_pivot);

	NNUnpairedModel(const std::string &_data_path, const PrimeStructure &_rna)
		: NNModel(_data_path, _rna) {}
	NNUnpairedModel(const std::string &_data_path)
		: NNModel(_data_path) {}

};

// TODO: Incorporate strain into the above model and have the algorithm figure it out.
/**
 * A special version of the NNUnpairedModel which incorporates strain.
 */
class StrainedUnpairedModel: public NNUnpairedModel {
	int strain = 30;
public:
	using NNUnpairedModel::MLClosure;
	/// The initiation cost of a multi-loop given only the relevant features.
	int MLClosure(int unpaired, int branches) const override;
	StrainedUnpairedModel(const std::string &_data_path, const PrimeStructure &_rna)
		: NNUnpairedModel(_data_path, _rna) {}
	StrainedUnpairedModel(const std::string &_data_path)
		: NNUnpairedModel(_data_path) {}
	librnary::energy_t MLStrain() const;
	void SetMLStrain(librnary::energy_t v);
};
}

#endif //RNARK_NN_UNPAIRED_MODEL_HPP

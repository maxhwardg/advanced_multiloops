//
// Created by max on 8/10/16.
//

#ifndef RNARK_AVERAGE_ASYM_MODEL_HPP
#define RNARK_AVERAGE_ASYM_MODEL_HPP

#include "nn_model.hpp"
#include "energy.hpp"

namespace librnary {
/**
 * Computes multi-loop energy using asymmetry in unpaired nucleotides to either side of a branch as a parameter.
 * By default, uses the parameters from the paper:
 * "Experimentally derived nearest-neighbor parameters for the stability of RNA three-and four-way multibranch loops"
 * by Mathews & Turner in 2002.
 */
class AverageAsymmetryModel: public NNModel {
	energy_t ml_init = 93, ml_branch = -6, ml_unpaired = 0, strain = 31;
	double ml_max_avg_asym = 2.0;
	kcalmol_t ml_asym_mult = 0.91;
public:
	int SumAsymmetry(const librnary::Surface &surf) const;
	energy_t MLClosure(const librnary::SSTree &sst, SSTreeNodeId id) const override;
	energy_t MLClosure(const librnary::Surface &surf) const override;
	energy_t MLClosure(int branches, int unpaired, int sum_asymmetry) const;
	energy_t MLClosureAsymCost(int sum_asymmetry, int branches) const;
	energy_t Strain() const;
	energy_t MLInit() const;
	energy_t MLBranchCost() const;
	energy_t MLUnpairedCost() const;
	kcalmol_t MLAsymmetryCoeff() const;
	double MLMaxAvgAsymmetry() const;
	void SetMLInit(energy_t v);
	void SetMLBranchCost(energy_t v);
	void SetMLUnpairedCost(energy_t v);
	void SetMLMaxAvgAsymmetry(double v);
	void SetMLAsymmetryCoeff(kcalmol_t v);
	void SetStrain(energy_t v);
	void SetMLParams(energy_t ml_init,
					 energy_t ml_branch,
					 energy_t ml_unpaired,
					 energy_t strain,
					 double ml_max_avg_asym,
					 kcalmol_t ml_asym_mult);
	AverageAsymmetryModel(const std::string &data_path, const PrimeStructure &_rna)
		: NNModel(data_path, _rna) {}
	AverageAsymmetryModel(const std::string &data_path)
		: NNModel(data_path) {}
};
}

#endif //RNARK_AVERAGE_ASYM_MODEL_HPP

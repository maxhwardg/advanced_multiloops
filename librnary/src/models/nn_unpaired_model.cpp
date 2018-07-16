//
// Created by max on 6/19/16.
//

#include "models/nn_unpaired_model.hpp"

librnary::energy_t librnary::NNUnpairedModel::MLBranchCost() const {
	return ml_br_cost;
}

librnary::energy_t librnary::NNUnpairedModel::MLClosure(int unpaired) const {
	return MLClosure(unpaired, 0);
}
librnary::energy_t librnary::NNUnpairedModel::MLClosure(int unpaired, int branches) const {
	assert(branches >= 0);
	assert(unpaired >= 0);
	int branch_cost = branches * MLBranchCost();
	if (unpaired <= ml_up_pivot)
		return ml_init + branch_cost + unpaired * ml_up_cost;
	else // TODO D*(int)round(log(...)) to get same behaviour as RNAstructure.
		return ml_init + branch_cost + ml_up_pivot * ml_up_cost
			+ librnary::KCalToEnergy(ml_log_mult * std::log(unpaired / (double) ml_up_pivot));
}

librnary::energy_t librnary::NNUnpairedModel::MLClosure(const librnary::SSTree &sstree,
														librnary::SSTreeNodeId node_id) const {
	int unpaired = sstree.Unpaired(node_id), branches = sstree.NumChildren(node_id) + 1;
	return MLClosure(unpaired, branches);

}
librnary::energy_t librnary::NNUnpairedModel::MLClosure(const librnary::Surface &surf) const {
	int unpaired = surf.Unpaired(), branches = surf.NumChildren() + 1;
	return MLClosure(unpaired, branches);
}

void librnary::NNUnpairedModel::SetMLParams(librnary::energy_t _ml_init,
											librnary::energy_t _ml_br_cost,
											librnary::energy_t _ml_up_cost,
											librnary::kcalmol_t _ml_log_mult,
											int _ml_up_pivot) {
	this->ml_init = _ml_init;
	this->ml_br_cost = _ml_br_cost;
	this->ml_up_cost = _ml_up_cost;
	this->ml_log_mult = _ml_log_mult;
	this->ml_up_pivot = _ml_up_pivot;
}

librnary::energy_t librnary::NNUnpairedModel::MLInitConstant() const {
	return ml_init;
}
void librnary::NNUnpairedModel::SetMLInitConstant(librnary::energy_t v) {
	NNUnpairedModel::ml_init = v;
}
void librnary::NNUnpairedModel::SetMLBranchCost(librnary::energy_t v) {
	NNUnpairedModel::ml_br_cost = v;
}
librnary::energy_t librnary::NNUnpairedModel::MLUnpairedCost() const {
	return ml_up_cost;
}
void librnary::NNUnpairedModel::SetMLUnpairedCost(librnary::energy_t v) {
	NNUnpairedModel::ml_up_cost = v;
}
librnary::kcalmol_t librnary::NNUnpairedModel::MLLogMultiplier() const {
	return ml_log_mult;
}
void librnary::NNUnpairedModel::SetMLLogMultiplier(librnary::kcalmol_t v) {
	NNUnpairedModel::ml_log_mult = v;
}
int librnary::NNUnpairedModel::MLUnpairedPivot() const {
	return ml_up_pivot;
}
void librnary::NNUnpairedModel::SetMLUnpairedPivot(int v) {
	NNUnpairedModel::ml_up_pivot = v;
}

librnary::energy_t librnary::StrainedUnpairedModel::MLStrain() const {
	return strain;
}

void librnary::StrainedUnpairedModel::SetMLStrain(librnary::energy_t v) {
	strain = v;
}

librnary::energy_t librnary::StrainedUnpairedModel::MLClosure(int unpaired, int branches) const {
	int strain_cost = branches == 3 && unpaired < 2 ? strain : 0;
	return librnary::NNUnpairedModel::MLClosure(unpaired, branches) + strain_cost;
}


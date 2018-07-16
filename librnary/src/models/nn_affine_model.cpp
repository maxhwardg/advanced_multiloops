//
// Created by max on 6/15/16.
//

#include "models/nn_affine_model.hpp"

librnary::energy_t
librnary::NNAffineModel::MLClosure(const librnary::SSTree &sstree, librnary::SSTreeNodeId node_id) const {
	return MLClosure(sstree.GetSurface(node_id));
}

librnary::energy_t librnary::NNAffineModel::MLClosure(const librnary::Surface &surf) const {
	return MLClosure(static_cast<unsigned>(surf.NumChildren() + 1), static_cast<unsigned>(surf.Unpaired()));
}

librnary::energy_t librnary::NNAffineModel::MLClosure(int branches, int unpaired) const {
	assert(branches >= 3);
	assert(unpaired >= 0);
	return MLInitCost() + MLUnpairedCost() * unpaired + MLBranchCost() * branches;
}


librnary::energy_t librnary::NNAffineModel::MLInitCost() const {
	return ml_init;
}

librnary::energy_t librnary::NNAffineModel::MLBranchCost() const {
	return ml_branch;
}

librnary::energy_t librnary::NNAffineModel::MLUnpairedCost() const {
	return ml_unpiared;
}

void librnary::NNAffineModel::SetMLInitCost(librnary::energy_t v) {
	ml_init = v;
}
void librnary::NNAffineModel::SetMLBranchCost(librnary::energy_t v) {
	ml_branch = v;
}
void librnary::NNAffineModel::SetMLUnpairedCost(librnary::energy_t v) {
	ml_unpiared = v;
}
void
librnary::NNAffineModel::SetMLParams(librnary::energy_t init, librnary::energy_t branch, librnary::energy_t unpaired) {
	SetMLInitCost(init);
	SetMLBranchCost(branch);
	SetMLUnpairedCost(unpaired);
}
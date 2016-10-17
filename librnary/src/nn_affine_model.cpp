//
// Created by max on 6/15/16.
//

#include "nn_affine_model.hpp"

librnary::energy_t
librnary::NNAffineModel::MLClosure(const librnary::SSTree &sstree, librnary::SSTreeNodeId node_id) const {
	return MLClosure(sstree.GetSurface(node_id));
}
librnary::energy_t librnary::NNAffineModel::MLClosure(const librnary::Surface &surf) const {
	return MLInitCost() + MLUnpairedCost() * surf.Unpaired() + MLBranchCost() * (surf.NumChildren() + 1);
}


librnary::energy_t librnary::NNAffineModel::MLInitCost() const {
	return dt->efn2a;
}

librnary::energy_t librnary::NNAffineModel::MLBranchCost() const {
	return dt->efn2c;
}

librnary::energy_t librnary::NNAffineModel::MLUnpairedCost() const {
	return dt->efn2b;
}
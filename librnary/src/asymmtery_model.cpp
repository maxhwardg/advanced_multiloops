//
// Created by max on 6/07/16.
//

#include "asymmetry_model.hpp"

librnary::energy_t librnary::AsymmetryModel::MLInit() const {
	return ml_init;
}
librnary::energy_t librnary::AsymmetryModel::MLBranchCost() const {
	return ml_branch;
}
librnary::energy_t librnary::AsymmetryModel::MLUnpairedCost() const {
	return ml_unpaired;
}
librnary::energy_t librnary::AsymmetryModel::MLAsymmetryCost() const {
	return ml_asymmetry;
}
void librnary::AsymmetryModel::SetMLInit(librnary::energy_t v) {
	ml_init = v;
}
void librnary::AsymmetryModel::SetMLBranchCost(librnary::energy_t v) {
	ml_branch = v;
}
void librnary::AsymmetryModel::SetMLUnpairedCost(librnary::energy_t v) {
	ml_unpaired = v;
}
void librnary::AsymmetryModel::SetMLAsymmetryCost(librnary::energy_t v) {
	ml_asymmetry = v;
}
void librnary::AsymmetryModel::SetMLParams(librnary::energy_t _ml_init,
										   librnary::energy_t _ml_branch,
										   librnary::energy_t _ml_unpaired,
										   librnary::energy_t _ml_asymmetry) {
	this->ml_init = _ml_init;
	this->ml_branch = _ml_branch;
	this->ml_unpaired = _ml_unpaired;
	this->ml_asymmetry = _ml_asymmetry;
}

int librnary::AsymmetryModel::SumAsymmetry(const librnary::Surface &surf) const {
	int sum_asymmetry = 0;
	int first_up = surf.Child(0).PairI() - surf.PairI() - 1;
	int prev_up = first_up;
	for (int i = 0; i < surf.NumChildren() - 1; ++i) {
		int next_up = surf.Child(i + 1).PairI() - surf.Child(i).PairJ() - 1;
		sum_asymmetry += abs(prev_up - next_up);
		prev_up = next_up;
	}
	int next_up = surf.PairJ() - surf.Child(surf.NumChildren() - 1).PairJ() - 1;
	sum_asymmetry += abs(prev_up - next_up);
	sum_asymmetry += abs(next_up - first_up);
	return sum_asymmetry;
}

librnary::energy_t librnary::AsymmetryModel::MLClosure(const librnary::SSTree &sst, SSTreeNodeId id) const {
	return MLClosure(sst.GetSurface(id));
}

librnary::energy_t librnary::AsymmetryModel::MLClosure(const librnary::Surface &surf) const {
	return MLClosure(surf.NumChildren() + 1, surf.Unpaired(), SumAsymmetry(surf));

}

librnary::energy_t librnary::AsymmetryModel::MLClosure(librnary::energy_t branches,
													   librnary::energy_t unpaired,
													   librnary::energy_t asymmetry) const {
	return ml_init + branches * ml_branch + unpaired * ml_unpaired + asymmetry * ml_asymmetry;
}
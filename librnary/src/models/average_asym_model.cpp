//
// Created by max on 8/11/16.
//

#include "models/average_asym_model.hpp"


using namespace std;

librnary::energy_t librnary::AverageAsymmetryModel::MLInit() const {
	return ml_init;
}
librnary::energy_t librnary::AverageAsymmetryModel::MLBranchCost() const {
	return ml_branch;
}
librnary::energy_t librnary::AverageAsymmetryModel::MLUnpairedCost() const {
	return ml_unpaired;
}

librnary::energy_t librnary::AverageAsymmetryModel::Strain() const {
	return strain;
}
void librnary::AverageAsymmetryModel::SetMLInit(int v) {
	ml_init = v;
}
void librnary::AverageAsymmetryModel::SetMLBranchCost(int v) {
	ml_branch = v;
}
void librnary::AverageAsymmetryModel::SetMLUnpairedCost(int v) {
	ml_unpaired = v;
}

double librnary::AverageAsymmetryModel::MLAsymmetryCoeff() const {
	return ml_asym_mult;
}
double librnary::AverageAsymmetryModel::MLMaxAvgAsymmetry() const {
	return ml_max_avg_asym;
}

void librnary::AverageAsymmetryModel::SetMLMaxAvgAsymmetry(double v) {
	ml_max_avg_asym = v;
}
void librnary::AverageAsymmetryModel::SetMLAsymmetryCoeff(kcalmol_t v) {
	ml_asym_mult = v;
}

void librnary::AverageAsymmetryModel::SetStrain(librnary::energy_t v) {
	strain = v;
}

void librnary::AverageAsymmetryModel::SetMLParams(energy_t _ml_init,
												  energy_t _ml_branch,
												  energy_t _ml_unpaired,
												  energy_t _strain,
												  double _ml_max_avg_asym,
												  kcalmol_t _ml_asym_mult) {
	this->ml_init = _ml_init;
	this->ml_branch = _ml_branch;
	this->ml_unpaired = _ml_unpaired;
	this->strain = _strain;
	this->ml_max_avg_asym = _ml_max_avg_asym;
	this->ml_asym_mult = _ml_asym_mult;
}

int librnary::AverageAsymmetryModel::SumAsymmetry(const librnary::Surface &surf) const {
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

librnary::energy_t librnary::AverageAsymmetryModel::MLClosure(const librnary::SSTree &sst, SSTreeNodeId id) const {
	return MLClosure(sst.GetSurface(id));
}

librnary::energy_t librnary::AverageAsymmetryModel::MLClosure(const librnary::Surface &surf) const {
	return MLClosure(surf.NumChildren() + 1, surf.Unpaired(), SumAsymmetry(surf));

}

librnary::energy_t librnary::AverageAsymmetryModel::MLClosure(int branches, int unpaired, int sum_asymmetry) const {
	return ml_init + branches * ml_branch + unpaired * ml_unpaired + MLClosureAsymCost(sum_asymmetry, branches)
		+ (branches == 3 && unpaired < 2 ? strain : 0);
}

librnary::energy_t librnary::AverageAsymmetryModel::MLClosureAsymCost(int sum_asymmetry, int branches) const {
	return KCalToEnergy(ml_asym_mult * min(ml_max_avg_asym, sum_asymmetry / static_cast<double>(branches)));
}
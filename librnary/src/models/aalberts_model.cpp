//
// Created by max on 6/22/16.
//

#include "models/aalberts_model.hpp"

using namespace std;

librnary::energy_t librnary::AalbertsModel::MLInit(int N, int M) const {
	// log_mult*ln(N^(6/5)*a^2 + M^(6/5)*b^2) + C
	// Multiply by ten to get tenths of kcal/mol.
	// Be careful to round AFTER this.
	assert(pow(N, power) * a * a + pow(M, power) * b * b >= 0);
	double inner_log = pow(N, power) * a * a + pow(M, power) * b * b;
	if (inner_log <= 0)
		return KCalToEnergy(C);
	else
		return KCalToEnergy(log_mult * log(pow(N, power) * a * a + pow(M, power) * b * b) + C);
}

// The initiation cost of a multi-loop given only the relevant features.
librnary::energy_t librnary::AalbertsModel::MLInitUpBr(int unpaired, int branches) const {
	int M = branches;
	int N = unpaired + M;
	return MLInit(N, M);
}

// The multi-loop initiation cost of a given loop region. Assumed to be a multi-loop.
librnary::energy_t
librnary::AalbertsModel::MLClosure(const librnary::SSTree &sstree, librnary::SSTreeNodeId node_id) const {
	int unpaired = sstree.Unpaired(node_id), branches = sstree.NumChildren(node_id) + 1;
	return MLInitUpBr(unpaired, branches);
}
// The multi-loop initiation cost of a given loop region. Assumed to be a multi-loop.
librnary::energy_t librnary::AalbertsModel::MLClosure(const librnary::Surface &surf) const {
	int unpaired = surf.Unpaired(), branches = surf.NumChildren() + 1;
	return MLInitUpBr(unpaired, branches);
}
void librnary::AalbertsModel::SetLogMultiplier(librnary::kcalmol_t v) {
	log_mult = v;
}
librnary::kcalmol_t librnary::AalbertsModel::LogMultiplier() const {
	return log_mult;
}
void librnary::AalbertsModel::SetAdditiveConstant(librnary::kcalmol_t v) {
	C = v;
}
librnary::kcalmol_t librnary::AalbertsModel::AdditiveConstant() const {
	return C;
}
void librnary::AalbertsModel::SetNCoeffBase(librnary::kcalmol_t v) {
	a = v;
}
double librnary::AalbertsModel::NCoeffBase() const {
	return a;
}
void librnary::AalbertsModel::SetMCoeffBase(double v) {
	b = v;
}
double librnary::AalbertsModel::MCoeffBase() const {
	return b;
}
void librnary::AalbertsModel::SetPower(double v) {
	power = v;
}
double librnary::AalbertsModel::Power() const {
	return power;
}
void librnary::AalbertsModel::SetMLParams(librnary::kcalmol_t _log_mult,
										  librnary::kcalmol_t _C,
										  double _a,
										  double _b,
										  double _power) {
	this->log_mult = _log_mult;
	this->C = _C;
	this->a = _a;
	this->b = _b;
	this->power = _power;
}
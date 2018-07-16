//
// Created by max on 6/22/16.
//


#ifndef RNARK_AALBERTS_MODEL_HPP
#define RNARK_AALBERTS_MODEL_HPP

#include "nn_model.hpp"
namespace librnary {
/**
 * This class implements the Aalbert & Nandagopal (2010) model of multi-loop energy.
 */
class AalbertsModel: public NNModel {
	/// The parameters used in energy calculation. Defaults that that found in the original publication.
	librnary::kcalmol_t log_mult = (59 / 36.0) * 0.61633135161, C = 0.0;
	double a = 6.2, b = 15, power = 6.0 / 5;
public:
	/// The initiation cost of a multi-loop given only the relevant features.
	energy_t MLInit(int N, int M) const;
	/// The initiation cost of a multi-loop given only the relevant features.
	energy_t MLInitUpBr(int unpaired, int branches) const;
	/// The multi-loop initiation cost of a given loop region. Assumed to be a multi-loop.
	energy_t MLClosure(const librnary::SSTree &sstree, librnary::SSTreeNodeId node_id) const override;
	/// The multi-loop initiation cost of a given loop region. Assumed to be a multi-loop.
	energy_t MLClosure(const librnary::Surface &surf) const override;
	void SetLogMultiplier(double v);
	librnary::kcalmol_t LogMultiplier() const;
	void SetAdditiveConstant(double v);
	librnary::kcalmol_t AdditiveConstant() const;
	void SetNCoeffBase(double v);
	double NCoeffBase() const;
	void SetMCoeffBase(double v);
	double MCoeffBase() const;
	void SetPower(double v);
	double Power() const;
	void SetMLParams(librnary::kcalmol_t log_mult,
					 librnary::kcalmol_t C,
					 double a,
					 double b,
					 double power);
	AalbertsModel(const std::string &data_path, const PrimeStructure &_rna)
		: NNModel(data_path, _rna) {}
	AalbertsModel(const std::string &data_path)
		: NNModel(data_path) {}
};
}

#endif //RNARK_AALBERTS_MODEL_HPP

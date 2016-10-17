//
// Created by max on 18/08/16.
//

#ifndef RNARK_AVERAGE_ASYM_SCORER_HPP
#define RNARK_AVERAGE_ASYM_SCORER_HPP

#include "nn_scorer.hpp"
#include "average_asym_model.hpp"


namespace librnary {
class AverageAsymmetryScorer: public NNScorer<AverageAsymmetryModel> {

public:
	virtual std::tuple<energy_t, librnary::energy_t> OptimalMLConfig(const Surface &surf) const override {
		using namespace std;

		auto ml_conf_vanilla = NNScorer<AverageAsymmetryModel>::OptimalMLConfig(surf);
		return make_tuple(get<0>(ml_conf_vanilla), em.MLClosure(surf));
	}
	AverageAsymmetryScorer(AverageAsymmetryModel _em) : NNScorer<AverageAsymmetryModel>(_em) {}

	AverageAsymmetryScorer(AverageAsymmetryModel _em, bool _stacking)
		: NNScorer<AverageAsymmetryModel>(_em, _stacking) {}
};
}

#endif //RNARK_AVERAGE_ASYM_SCORER_HPP

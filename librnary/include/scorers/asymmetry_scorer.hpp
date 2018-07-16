//
// Created by max on 7/7/16.
//

#include "nn_scorer.hpp"
#include "models/asymmetry_model.hpp"

#ifndef RNARK_ASYMMETRY_SCORER_HPP
#define RNARK_ASYMMETRY_SCORER_HPP

namespace librnary {
class AsymmetryScorer: public NNScorer<AsymmetryModel> {

public:
	virtual std::tuple<energy_t, librnary::energy_t> OptimalMLConfig(const Surface &surf) const override {
		using namespace std;

		auto ml_conf_vanilla = NNScorer<AsymmetryModel>::OptimalMLConfig(surf);
		return make_tuple(get<0>(ml_conf_vanilla), em.MLClosure(surf));
	}
	AsymmetryScorer(AsymmetryModel _em) : NNScorer<AsymmetryModel>(_em) {}

	AsymmetryScorer(AsymmetryModel _em, bool _stacking)
		: NNScorer<AsymmetryModel>(_em, _stacking) {}
};
}

#endif //RNARK_ASYMMETRY_SCORER_HPP

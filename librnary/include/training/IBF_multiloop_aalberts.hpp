//
// Created by max on 8/30/17.
//

#ifndef RNARK_IBF_MULTILOOP_AALBERTS_HPP
#define RNARK_IBF_MULTILOOP_AALBERTS_HPP
#include <models/aalberts_model.hpp>
#include <scorers/aalberts_scorer.hpp>

#include "IBF_multiloop.hpp"
namespace librnary {
/**
 * The Aalberts model is a special snowflake because stacking is affected by multi-loop parameters.
 * As such, it needs a weird ParamSetT interface, thus this class.
 */
template<typename ParamSetT>
class IBFMultiLoopAalberts: public IBFMultiLoop<ParamSetT, AalbertsModel, AalbertsScorer, AalbertsFolder> {

protected:

	/**
	 * Processes whatever is in fold results.
	 * This involves storing the structural information and energy.
	 * @return The average f-score of the fold results.
	 */
	double ProcessFoldResults() override {
		double sum_f_score_avgs = 0;
		for (size_t ctg = 0; ctg < this->cts.size(); ++ctg) {
			double sum_f_scores = 0;
			for (size_t i = 0; i < this->cts[ctg].size(); ++i) {
				double fscore = librnary::F1Score(this->fold_results[ctg][i], this->cts[ctg][i].match);
				sum_f_scores += fscore;
				if (this->fold_results[ctg][i] != this->cts[ctg][i].match
					&& this->false_multi_sets[ctg][i].count(this->fold_results[ctg][i]) == 0) {
					this->false_fscores[ctg][i].push_back(fscore);
					this->false_multi_sets[ctg][i].insert(this->fold_results[ctg][i]);
					librnary::SSTree sst(this->fold_results[ctg][i]);
					this->zero_ml_scorer.SetRNA(this->cts[ctg][i].primary);
					this->false_multi_base_energies[ctg][i].push_back(
						this->zero_ml_scorer.ScoreExterior(sst.RootSurface()));

					this->false_multi_info[ctg][i].emplace_back();

					std::vector<librnary::Surface> loops;
					librnary::ExtractMultiLoopSurfaces(loops, sst.RootSurface());
					for (const auto &loop : loops) {
						this->false_multi_info[ctg][i].back().emplace_back(this->zero_ml_scorer, loop);
					}
				}
			}
			sum_f_score_avgs += sum_f_scores / this->cts[ctg].size();
		}
		return sum_f_score_avgs / this->cts.size();
	}


public:

	/**
	 * @param _zero_model An energy model parameterization that gives multi-loops zero FE.
	 * @param _cts The list of CTs to use for training.
	 * @param _log Stream to log output to.
	 */
	IBFMultiLoopAalberts(const AalbertsModel &_zero_model, VV<CTData> _cts, std::ostream &_log)
		: IBFMultiLoop<ParamSetT, AalbertsModel, AalbertsScorer, AalbertsFolder>(_zero_model, _log) {
		using namespace std;
		this->cts = std::move(_cts);
		this->true_multi_info.resize(this->cts.size());
		this->true_multi_base_energies.resize(this->cts.size());
		for (size_t ctg = 0; ctg < this->cts.size(); ++ctg) {
			this->true_multi_base_energies[ctg].resize(this->cts[ctg].size());
			this->true_multi_info[ctg].resize(this->cts[ctg].size());
			for (size_t i = 0; i < this->cts[ctg].size(); ++i) {
				// Remove pseudoknots because they break everything.
				this->cts[ctg][i].match = librnary::RemovePseudoknotsMaximizePairs(this->cts[ctg][i].match);
				SSTree sst(this->cts[ctg][i].match);
				vector<librnary::Surface> loops;
				librnary::ExtractMultiLoopSurfaces(loops, sst.RootSurface());
				this->zero_ml_scorer.SetRNA(this->cts[ctg][i].primary);
				for (const auto &loop : loops) {
					this->true_multi_info[ctg][i].emplace_back(this->zero_ml_scorer, loop);
				}
				// Save scores without multi-loop score.
				this->true_multi_base_energies[ctg][i] = this->zero_ml_scorer.ScoreExterior(sst.RootSurface());
			}
		}
	}
};

}

#endif //RNARK_IBF_MULTILOOP_AALBERTS_HPP

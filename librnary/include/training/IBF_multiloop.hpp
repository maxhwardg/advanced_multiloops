//
// Created by max on 7/27/17.
//
// Contains an implementation of the iterative brute force method of training multi-loop parameters.
// This is similar to the iterative constraint generation method of Andronescu et al., but uses brute force instead of
// constraint programming for optimization. Only works well for small parameter spaces.

#ifndef RNARK_IBFML_HPP
#define RNARK_IBFML_HPP

#include <vector>
#include <iostream>
#include <set>

#include "read_cts.hpp"
#include "energy.hpp"
#include "ss_tree.hpp"
#include "multi_loop.hpp"
#include "statistics.hpp"
#include "parallel.hpp"
#include "vector_types.hpp"
#include "pseudoknot_removal.hpp"

namespace librnary {

/**
 * Iterative brute force method for training multi-loop parameters. Uses predicted F-score as a guide to fitness.
 * @tparam ParamSetT Parameter set type.
 * @tparam ModelT Energy model type.
 * @tparam ScorerT Energy model socrer type.
 */
template<typename ParamSetT, typename ModelT, typename ScorerT, typename FolderT>
class IBFMultiLoop {
protected:
	ModelT zero_model;
	VV<CTData> cts;
	std::ostream &log_stream;
	ScorerT zero_ml_scorer;
	VVV<typename ParamSetT::MultiInfo> true_multi_info;
	VVE true_multi_base_energies;

	VVV<double> false_fscores;
	_4DV<typename ParamSetT::MultiInfo> false_multi_info;
	VVVE false_multi_base_energies;

	std::vector<double> param_scores;

	VV<std::set<librnary::Matching>> false_multi_sets;
	VV<Matching> fold_results;

	size_t threads = std::thread::hardware_concurrency();

	int num_seeds = 0;
	int random_seed = 0;

	IBFMultiLoop(const ModelT &_zero_model, std::ostream &stream)
		: zero_model(_zero_model), log_stream(stream), zero_ml_scorer(zero_model) {}

	virtual long FindBestParams(const std::vector<ParamSetT> &params) {
		librnary::parallel_transform(params, param_scores, [=](const ParamSetT &pset) {
			auto local_model = zero_model;
			pset.LoadInto(local_model);

			double sum_averages = 0;

			for (size_t ctg = 0; ctg < cts.size(); ++ctg) {
				VE true_energies = true_multi_base_energies[ctg];
				auto min_base_energies = true_energies;

				V<int> min_energy_choice(cts[ctg].size(), -1);

				for (size_t i = 0; i < cts[ctg].size(); ++i) {
					for (const auto &mi : true_multi_info[ctg][i]) {
						true_energies[i] += mi.MLClosure(local_model);
					}
				}

				std::vector<energy_t> min_energies = true_energies;
				for (size_t i = 0; i < cts[ctg].size(); ++i) {
					for (size_t j = 0; j < false_multi_info[ctg][i].size(); ++j) {
						energy_t e = false_multi_base_energies[ctg][i][j];
						min_base_energies[i] = std::min(min_base_energies[i], e);
						for (const auto &mi : false_multi_info[ctg][i][j]) {
							e += mi.MLClosure(local_model);
						}
						// This is <= so that, if the true structure is MFE but non-unique, there is a penalty.
						if (e <= min_energies[i]) {
							min_energies[i] = e;
							min_energy_choice[i] = static_cast<int>(j);
						}
					}
				}
				double sum_fscores = 0;
				for (size_t i = 0; i < cts[ctg].size(); ++i) {
					if (min_energy_choice[i] == -1) {
						sum_fscores += 1;
					} else {
						sum_fscores += false_fscores[ctg][i][min_energy_choice[i]];
					}
				}
				sum_averages += sum_fscores / cts[ctg].size();
			}
			return sum_averages / cts.size();
		}, threads);

		auto it = librnary::parallel_max_element(begin(param_scores), end(param_scores), threads);
		return distance(begin(param_scores), it);
	}

	virtual void SeedStructures(const V<ParamSetT> &params, FolderT folder) {
		std::default_random_engine re(random_seed);
		for (int seed = 0; seed < num_seeds; ++seed) {
			auto param_set = params[re()%params.size()];
			FoldAllRNA(folder, param_set);
			double fscore = this->ProcessFoldResults();
			log_stream << "Seed #" << seed+1 << " " << param_set.to_string() << ": " << fscore << std::endl;
		}
	}

	virtual void InitTraining(const V<ParamSetT> &params, FolderT folder) {
		using namespace std;
		false_multi_sets = VV<set<librnary::Matching>>(cts.size());
		false_fscores = VVV<double>(cts.size());
		false_multi_info = _4DV<typename ParamSetT::MultiInfo>(cts.size());
		false_multi_base_energies = VVVE(cts.size());
		fold_results = VV<Matching>(cts.size());

		for (size_t ctg = 0; ctg < cts.size(); ++ctg) {
			false_multi_sets[ctg] = V<set<librnary::Matching>>(cts[ctg].size());
			false_fscores[ctg] = VV<double>(cts[ctg].size());
			false_multi_info[ctg] = VVV<typename ParamSetT::MultiInfo>(cts[ctg].size());
			false_multi_base_energies[ctg] = VVE(cts[ctg].size());
			fold_results[ctg] = V<Matching>(cts[ctg].size());
		}
		// Init arrays to be used repeatedly.
		param_scores = vector<double>(params.size());
		SeedStructures(params, folder);
	}

	virtual void FoldAllRNA(FolderT folder, ParamSetT param_set) {
		auto model = zero_model;
		param_set.LoadInto(model);
		for (size_t ctg = 0; ctg < cts.size(); ++ctg) {
			librnary::parallel_transform(cts[ctg], fold_results[ctg], [=](const librnary::CTData &ct_data) {
				auto local_folder = folder;
				local_folder.SetModel(model);
				local_folder.Fold(ct_data.primary);
				return local_folder.Traceback();
			}, threads);
		}
	}

	/**
	 * Processes whatever is in fold results.
	 * This involves storing the structural information and energy.
	 * @return The average f-score of the fold results.
	 */
	virtual double ProcessFoldResults() {
		double sum_f_score_avgs = 0;
		for (size_t ctg = 0; ctg < cts.size(); ++ctg) {
			double sum_f_scores = 0;
			for (size_t i = 0; i < cts[ctg].size(); ++i) {
				double fscore = librnary::F1Score(fold_results[ctg][i], cts[ctg][i].match);
				sum_f_scores += fscore;
				if (fold_results[ctg][i] != cts[ctg][i].match
					&& false_multi_sets[ctg][i].count(fold_results[ctg][i]) == 0) {
					false_fscores[ctg][i].push_back(fscore);
					false_multi_sets[ctg][i].insert(fold_results[ctg][i]);
					librnary::SSTree sst(fold_results[ctg][i]);
					zero_ml_scorer.SetRNA(cts[ctg][i].primary);
					false_multi_base_energies[ctg][i].push_back(zero_ml_scorer.ScoreExterior(sst.RootSurface()));

					false_multi_info[ctg][i].emplace_back();

					std::vector<librnary::Surface> loops;
					librnary::ExtractMultiLoopSurfaces(loops, sst.RootSurface());
					for (const auto &loop : loops) {
						false_multi_info[ctg][i].back().emplace_back(loop);
					}
				}
			}
			sum_f_score_avgs += sum_f_scores / cts[ctg].size();
		}
		return sum_f_score_avgs / cts.size();
	}


public:
	void SetNumStructureSeeds(int num) {
		assert(num >= 0);
		num_seeds = num;
	}
	void SetRandomSeed(int rnd) {
		random_seed = rnd;
	}
	void SetThreads(size_t num_threads) {
		threads = num_threads;
	}
	/**
	 * @param params List of parameters to optimize.
	 * @param init Starting parameter set.
	 * @param folder Folding algorithm used to make predictions.
	 * @tparam FolderT Folder type to use.
	 * @param num_epochs Number of training iterations until exit.
	 * @return The optimized parameter set.
	 */
	virtual ParamSetT Train(const V<ParamSetT> &params,
					ParamSetT init,
					FolderT folder,
					int num_epochs) {
		using namespace std;

		// Arbitrary parameter set to start with.
		ParamSetT param_set = init;

		double best_sc = -1;
		ParamSetT best_set = init;

		// Training loop.
		InitTraining(params, folder);


		for (int epoch = 0; epoch < num_epochs; ++epoch) {

			// Fold all the RNAs using the current parameter set.
			FoldAllRNA(folder, param_set);

			// Add new folds to the false_multi_loop lists (if they are incorrect).
			// Also compute the F-scores of the folds.
			double avg_f_score = this->ProcessFoldResults();

			if (avg_f_score > best_sc) {
				best_sc = avg_f_score;
				best_set = param_set;
			}

			log_stream << "Epoch #" << epoch << ": " << endl << "\tAverage F-Score = " << avg_f_score << endl;

			// Find a parameter set that minimises the RMSE over all fold results.

			long ind = FindBestParams(params);

			log_stream << "\tBest score = " << param_scores[ind] << endl;
			log_stream << "\t" << params[ind].to_string() << endl;


			// If we're stuck in a loop, break.
			if (params[ind] == param_set)
				break;
			param_set = params[ind];
		}
		return best_set;
	}
	/**
	 * @param _zero_model An energy model parameterization that gives multi-loops zero FE.
	 * @param _cts The list of CTs to use for training.
	 * @param _log Stream to log output to.
	 */
	IBFMultiLoop(const ModelT &_zero_model, VV<CTData> _cts, std::ostream &_log)
		: zero_model(_zero_model), cts(std::move(_cts)), log_stream(_log), zero_ml_scorer(zero_model) {
		using namespace std;
		true_multi_info.resize(cts.size());
		true_multi_base_energies.resize(cts.size());
		for (size_t ctg = 0; ctg < cts.size(); ++ctg) {
			true_multi_base_energies[ctg].resize(cts[ctg].size());
			true_multi_info[ctg].resize(cts[ctg].size());
			for (size_t i = 0; i < cts[ctg].size(); ++i) {
				// Remove pseudoknots because they break everything.
				cts[ctg][i].match = librnary::RemovePseudoknotsMaximizePairs(cts[ctg][i].match);
				SSTree sst(cts[ctg][i].match);
				vector<librnary::Surface> loops;
				librnary::ExtractMultiLoopSurfaces(loops, sst.RootSurface());
				for (const auto &loop : loops) {
					true_multi_info[ctg][i].emplace_back(loop);
				}
				// Save scores without multi-loop score.
				zero_ml_scorer.SetRNA(cts[ctg][i].primary);
				true_multi_base_energies[ctg][i] = zero_ml_scorer.ScoreExterior(sst.RootSurface());
			}
		}
	}
};
}

#endif //RNARK_IBFML_HPP

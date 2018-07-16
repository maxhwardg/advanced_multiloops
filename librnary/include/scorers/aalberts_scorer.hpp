//
// Created by max on 6/27/16.
//

#ifndef RNARK_AALBERTS_SCORER_HPP
#define RNARK_AALBERTS_SCORER_HPP

#include "nn_scorer.hpp"
#include "multi_array.hpp"
#include "models/aalberts_model.hpp"
namespace librnary {


class AalbertsScorer: public NNScorer<AalbertsModel> {
private:
	/**
	 * Returns a DP table containing the optimal configurations of stacking and a/b length segments.
	 * dp_table[as][bs][s] is the MFE configuration from subsurface s onward.
	 * The configuration contains as  a-length segments, and bs b-length segments, in the prefix preceding s.
	 * Optimal "stacking" is defined here http://rna.urmc.rochester.edu/NNDB/turner04/exterior.html
	 * Aalberts & Nandagopal (2010) define a and b length segments.
	 * Takes a list of features (branches/unpaired) as input. See MLFeatures(loop) for details.
	 */
	Array3D <energy_t> MakeStackingTable(const std::vector<StackingFeature> &features) const;

	/**
	 * Traces back multi-loop stacking.
	 * @param features Featues in the multi-loop.
	 * @param table Previously filled DP table (from MakeStackingTable).
	 * @param as Number of A segments to start with. Index into DP table.
	 * @param bs Number of B segments to start with. Index into DP table.
	 * @param ind Feature of multi-loop to start with. Index into DP table.
	 * @return A list of stacking interactions in the optimal stacking.
	 */
	std::vector<Stacking>
	TraceStackingTable(const std::vector<StackingFeature> &features,
					   const Array3D <energy_t> &table,
					   int as,
					   int bs,
					   int ind) const;

public:


	// TODO Override trace functions

	/**
	 * Calculates the FE of the optimal configuration of dangles for a multiloop surface.
	 * Also considers coaxial stacking and a/b length segments.
	 */
	std::tuple<energy_t, energy_t> OptimalMLConfig(const Surface &surf) const override;

	std::tuple<std::vector<Stacking>, ClosureFeaturesMap> TraceMLConfig(const Surface &surf) const override;




	explicit AalbertsScorer(const AalbertsModel &_em)
		: NNScorer(_em) {}

	AalbertsScorer(const AalbertsModel &_em, bool _stacking)
		: NNScorer(_em, _stacking) {}
};
}

#endif //RNARK_AALBERTS_SCORER_HPP

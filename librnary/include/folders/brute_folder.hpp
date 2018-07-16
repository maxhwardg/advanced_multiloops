//
// Created by max on 6/19/16.
//

#ifndef RNARK_BRUTE_FOLDER_HPP
#define RNARK_BRUTE_FOLDER_HPP

#include "ss_enumeration.hpp"
#include "ss_tree.hpp"

#include <queue>

namespace librnary {
/**
 * Finds MFE folds of an RNA sequence by brute force.
 */
class BruteFolder {
	StructureEnumerator enumer;
	/**
	 * FoldeMFEFunctor is used to ensure blazingly fast brute force searches.
	 * It avoids all virtual calls as would be done by using a lambda.
	 */
	template<typename EnergyModel>
	struct FoldMFEFunctor {
		EnergyModel em;
		librnary::energy_t mfe, delta;
		std::priority_queue<std::tuple<int, Matching>> mfe_structs;

		FoldMFEFunctor(const EnergyModel &_em, int _delta)
			: em(_em), delta(_delta) {
			mfe = em.MaxMFE();
		}
		void operator()(const Matching &m) {
			SSTree sstree(m);
			int sc = em.ScoreExterior(sstree.RootSurface());
			if (sc < mfe) {
				// New MFE; update mfe variable.
				mfe = sc;
				mfe_structs.push(make_tuple(sc, m));
				// Remove out of range structures. The queue is a max-first queue.
				while (std::get<0>(mfe_structs.top()) > mfe + delta)
					mfe_structs.pop();
			} else if (sc <= mfe + delta) // Push if within range.
				mfe_structs.push(make_tuple(sc, m));
		}

	};
public:
	explicit BruteFolder(const StructureEnumerator &_enumer)
		: enumer(_enumer) {}
	/**
	 * Given a primary sequence and the number of structures to find, returns a vector containing the nMFE lowest
	 * free energy structures.
	 * Each element in the vector is a tuple of free energy and the structure.
	 * These are found using an exhaustive search.
	 */
	template<typename EnergyModel>
	std::vector<std::tuple<librnary::energy_t, Matching>>
	FoldN(EnergyModel &em, const PrimeStructure &primary, unsigned nMFE) {
		em.SetRNA(primary);
		std::priority_queue<std::tuple<librnary::energy_t, Matching> > pq;
		auto f = [&pq, &primary, nMFE, &em](const Matching &m) {
			SSTree sstree(m);
			int sc = em.ScoreExterior(sstree.RootSurface());
			if (pq.size() < nMFE)
				pq.push(make_tuple(sc, m));
			else if (sc < std::get<0>(pq.top())) {
				pq.pop();
				pq.push(make_tuple(sc, m));
			}
		};
		enumer.Enumerate(primary, f);
		std::vector<std::tuple<librnary::energy_t, Matching>> best;
		while (!pq.empty()) {
			best.push_back(pq.top());
			pq.pop();
		}
		return best;
	}
	/**
	 * Takes a primary sequence. Returns the MFE score and a complete list of MFE structures.
	 * Also finds structures within delta of the MFE.
	 * Uses exhaustive search.
	 */
	template<typename EnergyModel>
	std::vector<std::tuple<librnary::energy_t, Matching>>
	FoldMFE(EnergyModel &em, const PrimeStructure &primary, librnary::energy_t delta = 0) {
		em.SetRNA(primary);
		FoldMFEFunctor<EnergyModel> func(em, delta);
		enumer.Enumerate(primary, func);
		std::vector<std::tuple<librnary::energy_t, Matching>> res;
		res.reserve(func.mfe_structs.size());
		while (!func.mfe_structs.empty()) {
			res.push_back(func.mfe_structs.top());
			func.mfe_structs.pop();
		}
		return res;
	}
};
}

#endif //RNARK_BRUTE_FOLDER_HPP

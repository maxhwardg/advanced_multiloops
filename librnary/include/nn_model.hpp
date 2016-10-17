//
// Created by max on 5/31/16.
//

/**
 * This file defines a class that implements that nearest neighbour energy model as described by Mathews & Turner,
 * Zuker et al., and Tinoco et al.
 */

#ifndef RNARK_NN_MODEL_HPP
#define RNARK_NN_MODEL_HPP

#include "rna_structure.hpp"
#include "ss_tree.hpp"
#include "energy.hpp"

namespace librnary {

/** An instance of NNModel represents a particular parametrization of the nearest neighbour energy model.
 * Currently based largely on the 2004 version found at http://rna.urmc.rochester.edu/NNDB/.
 * Because multi-loops are treated oddly (linear for DP, logarithmic for energy evaluation) they are abstracted out.
 * Note that all the energy functions are orthogonal, so, for example, no function will call Branch for you.
 * In addition, the functions provided have a philosophy of being as clever as possible. They will deduce what kind of
 * coaxial stack you mean--whether i,j closes a multi-loop--for example. */
class NNModel {
protected:
	PrimeStructure rna;
	std::shared_ptr<datatable> dt;
	std::shared_ptr<structure> struc;
public:
	/// Minimum number of unpaired nucleotides allowed in a hairpin loop.
	static const int MIN_HAIRPIN_UNPAIRED = 3;

	/** Returns the energy of a one-loop from 5' i to 3' j.
	 * A one-loop is a hairpin loop.
	 * Assumes i<j.
	 * http://rna.urmc.rochester.edu/NNDB/turner04/hairpin.html */
	virtual energy_t OneLoop(int i, int j) const;

	/** Returns the energy of a two-loop between 5' -> i -> k -> l -> j -> 3'.
	 * (_(_)_)
	 * i k l j
	 * This could be a bulge, stacking region, or internal loop.
	 * http://rna.urmc.rochester.edu/NNDB/turner04/bulge.html
	 * http://rna.urmc.rochester.edu/NNDB/turner04/wc.html (stacking)
	 * http://rna.urmc.rochester.edu/NNDB/turner04/internal.html */
	virtual energy_t TwoLoop(int i, int k, int l, int j) const;

	/** The free energy cost of a (external/multi)-loop branch closed by 5' i and 3' j.
	 * Considers AU/GU closure penalties.
	 * Note that this function doesn't care about i<j. */
	virtual energy_t Branch(int i, int j) const;

	/** Computes the initiation cost of a multi-loop. Returns zero since we consider multi-loops to have no model here.
	 * It is used to build other multiloop models polymorphically as derived classes.
	 * Note that this function assumes that the given node closes a valid multiloop. */
	virtual energy_t MLClosure(const librnary::SSTree &sstree, librnary::SSTreeNodeId node_id) const = 0;

	/** Computes the initiation cost of a multi-loop. Returns zero since we consider multi-loops to have no model here.
	 * It is used to build other multiloop models polymorphically as derived classes.
	 * Note that this function assumes that the given node closes a valid multiloop. */
	virtual energy_t MLClosure(const librnary::Surface &surf) const = 0;

	/** A flush coaxial stacking between pairs i,j and k,l.
	 * 5' i and 3' j are a bonding pair which are in a coaxial stack with 5' k and 3' l.
	 * Assumes i,j < k,l, or that i,j closes a multi-loop.
	 * Note that this only computes the coaxial stacking score.
	 * Thus Branch WILL NOT be called by this function.
	 * http://rna.urmc.rochester.edu/NNDB/turner04/coax.html */
	virtual energy_t FlushCoax(int i, int j, int k, int l) const;

	/** Mismatch mediated coaxial stack. Mismatched bases are next to stem at i,j.
	 * Parameters k and l represent the other stem in the coaxial stack.
	 * This function assumes that either i,j is fully left or right of k,l, or that i,j or k,l closes a multi-loop.
	 * Note that this only computes the coaxial stacking score.
	 * Thus Branch WILL NOT be called by this function.
	 * http://rna.urmc.rochester.edu/NNDB/turner04/coax.html */
	virtual energy_t MismatchCoax(int i, int j, int k, int l) const;

	/** Returns the free energy change of a 5' dangle off the pair i,j.
	 * Assumes i<j and that i,j doesn't close a multi-loop. */
	virtual energy_t FiveDangle(int i, int j) const;

	/**
	 * Same as FiveDangle, but assumes i,j close a multi-loop.
	 * As such, i<dangle<j.
	 */
	virtual energy_t ClosingFiveDangle(int i, int j) const;

	/** Returns the free energy change of a 3' dangle off the pair i,j.
	 * Assumes i<j and that i,j doesn't close a multi-loop. */
	virtual energy_t ThreeDangle(int i, int j) const;

	/**
	 * Same as ThreeDangle, but assumes i,j close a multi-loop.
	 * As such, i<dangle<j.
	 */
	virtual energy_t ClosingThreeDangle(int i, int j) const;

	/** Returns the free energy change of a terminal mismatch off the pair i,j.
	 * Assumes i<j and that i,j doesn't close a multi-loop. */
	virtual energy_t Mismatch(int i, int j) const;

	/**
	 * Same as Mismatch, but assumes i,j close a multi-loop.
	 * As such, i<mismatch<j.
	 */
	virtual energy_t ClosingMismatch(int i, int j) const;

	/// Returns a safe "infinity" value.
	virtual energy_t MaxMFE() const;

	virtual ~NNModel() {}

	void SetRNA(const PrimeStructure &rna);

	librnary::PrimeStructure RNA() const;

	NNModel(const std::string &data_path, const PrimeStructure &_rna) {
		this->dt = librnary::LoadDatatable(data_path);
		SetRNA(_rna);
	}

	NNModel(const std::string &data_path) {
		this->dt = librnary::LoadDatatable(data_path);
	}

};
}

#endif //RNARK_NN_MODEL_HPP

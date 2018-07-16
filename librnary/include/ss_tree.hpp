//
// Created by max on 5/5/16.
//

/**
 * This file contains the secondary structure tree class and related types and functions.
 */

#ifndef RNARK_SS_TREE_HPP
#define RNARK_SS_TREE_HPP

#include "secondary_structure.hpp"

namespace librnary {

/// Index type of a node in an SSTree. Should be able to operate like a positive integer.
typedef int SSTreeNodeId;

class SSTree;
/**
 * This represents a surface in the Rivas & Eddy (1999) sense.
 * Used to succinctly represent an arc and its subarcs. In other words, a node and its child edges in the SSTree.
 * This is the class that is designed for easy, const correct epxloration of an SSTree.
 * It should be used as an iterator for an SSTree. Be careful that the SSTree does not go out of scope!
 */
class Surface {
private:
	const SSTree *tree;
	SSTreeNodeId id;
public:
	Surface(const SSTree *_tree, SSTreeNodeId _id)
		: tree(_tree), id(_id) {}
	/// Returns the left/5' nucleotide in the arc closing the surface.
	int PairI() const;
	/// Returns the right/3' nucleotide in the arc closing the surface.
	int PairJ() const;
	/// Returns the number of unpaired nucleotides enclosed by the surface.
	int Unpaired() const;
	/// Returns the parent surface in the SSTree.
	Surface Parent() const;
	/// Create and return a list of children in the SSTree.
	std::vector<Surface> Children() const;
	/// Given the (zero indexed) index of a child, returns the corresponding surface.
	Surface Child(int idx) const;
	/// Returns the number of children.
	int NumChildren() const;
	/// Get the root surface containing this surface.
	Surface RootSurface() const;
	/// Returns true iff the surface is an external loop surface.
	bool IsExternalLoop() const;
	bool operator==(const Surface &rhs) const;
};

/**
 * Determines the loop type closed by a surface.
 * @param surf The surface.
 * @return The loop typ.e
 */
LoopType ClosedLoopType(const Surface &surf);


/**
 * Determines if a surface closes a stacking loop.
 * @param surf The surface to check.
 * @return True if surf closes a stacking region, and false otherwise.
 */
bool IsStacking(const Surface &surf);

/**
 * SSTree is a secondary structure tree rooted at a particular pairing.
 * Secondary structure can be represented as a tree.
 * This is because bonds are nested like parenthesise.
 * This class is essentially a thin wrapper around an adjacency list representation of a tree.
 * Ids are guaranteed to be dense [0, N) integer-like keys.
 *
 * Let us say that every arc in a pseudoknot free secondary structure is a node in a tree.
 * The children of these nodes are the accessible arcs in the Lyngso sense.
 * The root of the tree is a special arc encompassing all others.
 */
class SSTree {
	friend class Surface;
	/// Adjacency list representation of tree.
	std::vector<std::vector<SSTreeNodeId>> adj_list;
	/// Direct parent in the tree
	std::vector<SSTreeNodeId> parents;
	/// List of bonding pairs in the tree. Index corresponds to SSTreeNodeId.
	std::vector<BondPair> bond_pairs;
	SSTreeNodeId next_id;
	/// Range (i,j) is exclusive.
	void Construct(const Matching &matching, int i, int j);
	void Construct(const Matching &matching, SSTreeNodeId curr, int i, int j);
public:
	/// Construct from Matching.
	explicit SSTree(const Matching &matching);
	~SSTree() = default;
	/// Returns the number of unpaired nucleotides accessible from the arc-node.
	int Unpaired(SSTreeNodeId node_id) const;
	/// Whether a node-arc is the external loop.
	bool IsExternalLoop(SSTreeNodeId node_id) const;
	/// The number of children for a node.
	int NumChildren(SSTreeNodeId node_id) const;
	/// Returns the adjacency list for a particular node.
	std::vector<SSTreeNodeId> &Children(SSTreeNodeId node_id);
	/// A child's node id. Uses zero indexing and children will be in left to right order.
	SSTreeNodeId Child(SSTreeNodeId node_id, int child) const;
	int PairI(SSTreeNodeId node_id) const;
	int PairJ(SSTreeNodeId node_id) const;
	/// The id of the root arc-node.
	SSTreeNodeId RootId() const;
	/// Get a Surface representation of the root.
	Surface RootSurface() const;
	Surface GetSurface(SSTreeNodeId id) const;
	/// Returns the id of the parent arc-node in the tree.
	SSTreeNodeId Parent(SSTreeNodeId node_id) const;
	/// Returns the number of nodes. Also largest node id + 1.
	int NumNodes() const;
};

}

#endif //RNARK_SS_TREE_HPP

#pragma once
#ifndef SAFRA_TREE_H
#define SAFRA_TREE_H

#include <string>
#include <boost/dynamic_bitset/dynamic_bitset.hpp>
#include <boost/functional/hash.hpp>

#include "utils.hpp"
#include "NBW.hpp"


// SETTINGS FOR SAFRA'S CONSTRUCTION

/** Whether or not to mark newly created child nodes.
 */
#define MARK_NEW_CHILDREN true

/** Whether to transition the labels before creating children.
 */
#define TRANSITION_FIRST true

/* Whether to save the Safra trees in memory until the next 
 * determinization so that the data is still viewable (for example
 * with @function DRW::to_GASt_string() ).
 * If set to true, the trees themselves are only deleted (and the 
 * memory freed) with a call to @function SafraTree::reset().
 */
#define SAVE_TREE_DATA true

class SafraTree{
  private:
  
    /* Keeps track of all instances of SafraTrees. These instances are deleted
     * with a call to @function SafraTree::reset().
     */
    static std::vector<SafraTree*> references;


  public:
    static int next_tree_id;
    
    /** After a determinization, this vector contains references to the
     * trees (corresponding to the states of the automaton) if the value of
     * SAVE_TREE_DATA is true. Otherwise, it is empty. To clear it and free
     * any memory that is being held by determinization structures for
     * debugging purposes, call @function SafraTree::reset().
     */
    static std::vector<SafraTree*> trees;



    /* ------------- Fields ---------------- */
    
    int treeID; // globally unique

    int name; // the state that this represents in the Rabin automaton

    SafraNode* root;
    SafraTree** targets;
        
    // Precomputed hash value.
    std::size_t hvalue;

    boost::dynamic_bitset<unsigned long> marked_nodes;
    boost::dynamic_bitset<unsigned long> used_node_names;
    boost::dynamic_bitset<unsigned long> temp_node_names;
    
    SafraNode* node_storage; // a pointer to the only place that they reside in memory
  
    
    /* ------------- Methods --------------- */
    
    
    /*
     * Constructor: trying to do as little as possible here.
     * Doesn't do much, although you need to know the size of the NBW
     * we're determinizing
     */
    SafraTree(int buchi_size, int alphabet_size);
    ~SafraTree();    
    
    /** Reset and/or initialize static variables, before or between 
     *  determinizations.
     */
    static void reset();

    static SafraTree* build_initial_tree(const NBW& input_automaton);
    static SafraTree* get_transition(const SafraTree& old_tree, const NBW& input, int character);
    
    /**
     * Get the SafraTree corresponding to state @param i. Only works if
     * SAVE_TREE_DATA is true, since this accesses the private static vector
     * @field trees.
     */
    static SafraTree* get_tree(int i);

    bool operator==(const SafraTree& other) const;
    
    /** Get the lowest available node name for a new node and remove it
     * from the list.
     */
    int name_node();
    
    /** Return a node name to the list of available names -- do this when
     * destroying a node.
     */
    void free_node_name(int name);
    
    /** We'll maintain a list of marked nodes as well, to speed up calculating
     * the acceptance condition.
     */
    void mark_node(int name);
    
    /** Reserve a node name temporarily. Called when the child node would be
     * immediately deleted, but the name given to the child might influence
     * the naming of other nodes created before that deletion takes place.
     * This should be called with a name that was previously reserved by a 
     * call to @function name_node.
     */
    void temp_name_node(int name);
    
    /** Frees any node names which were "temporarily reserved" with 
     * a call to temp_name_node().
     */
    void free_temp_names();
    
    /** Generate a textual representation of the Safra tree.
     */
    std::string to_string() const;
    
    /** Because one to_string method isn't enough for some folks, this one
     * outputs the tree data in GASt format.
     */
    std::string to_GASt_string() const;
};

class SafraNode{

  private:
    /** Used only for differentiating instances, primarily in debugging.
     */
    int ID;
    /** For assigning IDs.
     */
    static int next_id;

  public:
    int name;
    state_set_t states;
    bool marked;
    
    /** These are the primary references to the other nodes in the tree, so
     * any deletion of nodes and freeing of memory recurses using this list.
     */
    std::vector<SafraNode*> children;
    
    /** A pointer to the tree is kept so that this node can tell the tree to
     * free its name when it is deleted. The node (obviously) does not uniquely
     * own the tree.
     */    
    SafraTree* tree;
    
    
    SafraNode();
    ~SafraNode();
    
    /** Frees the node name (and all of the children) recursively. Since memory
     * for the nodes is now stored in the tree no actual memory is freed.
     */
    void kill_node();
        		
    bool operator==(const SafraNode& other) const;
    
    /**
     * Updates a Safra tree using one depth-first search.
     * Copy nodes (erasing all marks), create children, and perform transitions of labels. 
     * Then suppress states, suppress node labelings, and mark appropriate nodes for
     * Safra's construction. 
     * @return the new root node of the tree, or NULL if the tree (or subtree) is empty.
     */
    SafraNode* clone_spawn_and_transition(SafraTree* new_tree, const NBW& input, int character, bool root, state_set_t& kill_set, const state_set_t& nbw_final_states);
    
    /** Generate a string representation of the node -- the indent level is used
     * since this function is recursive, so that the tree is printed with an
     * easy-to-interpret format.
     */
    std::string to_string(int indent_level) const;

    /** A slightly more verbose form of the to_string method which 
     * conforms more exactly to the GASt syntax.
     */
    std::string to_GASt_string(int indent_level) const;
    
    /** Generate a hash value for this node. Recursive.
     */
    std::size_t hash_value() const;
    
    /** Accumulate a list of names for all nodes in the subtree rooted at this node. Recursive.
     */
    void accumulate_subtree_names(state_set_t& subtree_names);
    
    friend class SafraTree;
};


/** A custom hash function for pointers to Safra trees.
 * Returns a value based on the hash of the underlying tree.
 */
struct stp_hash_t
    : std::unary_function<SafraTree*, std::size_t>
{
    std::size_t operator()(SafraTree* const& st_p) const
    {
        return st_p->hvalue;       
    }
};

/** A custom equality comparison for pointers to Safra trees.
 * Returns true if the underlying trees are equal.
 */
struct stp_eq_t
    : std::binary_function<SafraTree*, SafraTree*, bool>
{
    bool operator()(SafraTree* const& x,
        SafraTree* const& y) const
    {
        return (*x) == (*y);
    }
};


#endif
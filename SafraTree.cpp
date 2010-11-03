/* @file SafraTree.cpp
 * @author Joe Gershenson
 *
 * Implementation of functions for Safra trees, used in Safra's construction.
 * A SafraTree is composed of SafraNodes -- also defined in this class.
 */


#include <iostream>
#include <string>
#include <iterator>
#include <vector>
#include <boost/dynamic_bitset/dynamic_bitset.hpp>
#include <boost/functional/hash.hpp>

#include "SafraTree.hpp"

int SafraNode::next_id = 0;
int SafraTree::next_tree_id = 0;
std::vector<SafraTree*> SafraTree::trees;

std::vector<SafraTree*> SafraTree::references;

/*** Implementation of SafraTree ***/

SafraTree::SafraTree(int buchi_size, int alphabet_size){
    this->name = -1; //trees start unnamed.
    this->treeID = next_tree_id++;
    
    this->targets = new SafraTree*[alphabet_size];
    
    this->marked_nodes.resize(2*buchi_size);
    this->used_node_names.resize(2*buchi_size);
    this->temp_node_names.resize(2*buchi_size);
    
    this->node_storage = new SafraNode[2*buchi_size];
    for (int i = 0; i < 2*buchi_size; i++)
        node_storage[i].tree = NULL;
        
    references.push_back(this);
}

SafraTree::~SafraTree(){
    delete[] this->targets;
    delete[] this->node_storage;
}

void SafraTree::reset(){    
    SafraTree::next_tree_id = 0;
    SafraNode::next_id = 0;
    // delete old references to free memory 
    for(int i = 0; i < SafraTree::references.size(); i++)
        delete SafraTree::references[i];
    SafraTree::references.clear();
    SafraTree::trees.clear();
}

bool SafraTree::operator==(const SafraTree& other) const {
    if((this->root == NULL) || (other.root == NULL))
        return (this->root == NULL) && (other.root == NULL);
    else{
        if(this->hvalue != other.hvalue)
            return false;
        if(this->used_node_names != other.used_node_names)
            return false;
        return (*(this->root) == *(other.root));        
    }
}

int SafraTree::name_node(){
    int i = 0;
    while(this->used_node_names[i] == true){
        i++;
    }     
    this->used_node_names.set(i);
    return (i+1);
}

void SafraTree::free_node_name(int name){
    this->used_node_names.reset(name-1);
    this->marked_nodes.reset(name-1);
}

inline void SafraTree::mark_node(int name){
    this->marked_nodes.set(name-1);
}

inline void SafraTree::temp_name_node(int name){ 
    this->temp_node_names.set(name-1);
}

inline void SafraTree::free_temp_names(){
    this->used_node_names -= this->temp_node_names;
    this->temp_node_names.reset();
}

SafraTree* SafraTree::build_initial_tree(const NBW& input){
    state_set_t nbw_initial_states = input.get_initial_states();
    state_set_t nbw_final_states = input.get_final_states();

    SafraTree* ret = new SafraTree(input.size, input.alphabet_size);

    ret->name = 0;

    // construct new root node
    ret->root = new(ret->node_storage) SafraNode();
    ret->root->ID = SafraNode::next_id++;
    ret->root->name = ret->name_node();
    ret->root->tree = ret;
    ret->root->states = nbw_initial_states;
    
    state_set_t copy(nbw_initial_states);
    copy &= nbw_final_states;
    if(copy.none()){ // true if no bits are set, so the sets are disjoint
        ret->root->marked = false;
    } else if (nbw_initial_states.is_subset_of(nbw_final_states)) {
        ret->root->marked = true;
    } else {
        ret->root->marked = false;
        
        void* new_node_loc = (void*)(ret->node_storage + 1);
        
        SafraNode* child = new(new_node_loc) SafraNode();
        child->ID = SafraNode::next_id++;
        child->name = ret->name_node();
        child->tree = ret;
        child->states = copy;
        child->marked = true;
        ret->root->children.push_back(child);
    }
    
    // precompute hash value
    ret->hvalue = ret->root->hash_value();
    
    return ret;    
}

/** 
 * Clone and transition a Safra tree.
 *
 */
SafraTree* SafraTree::get_transition(const SafraTree& old_tree, const NBW& input, int character){
    if(old_tree.root == NULL){
        SafraTree* ret = new SafraTree(input.size, input.alphabet_size);
        ret->root = NULL;
        ret->hvalue = 0;
        return ret;
    } else {
        SafraTree* ret = new SafraTree(input.size, input.alphabet_size);
        ret->used_node_names = old_tree.used_node_names;
        
        state_set_t kill_set(input.size);
               
        old_tree.root->clone_spawn_and_transition(ret, input, character, true, kill_set, input.get_final_states());

        // precompute hash value
        if(ret->root == NULL){
            ret->hvalue = 0;
        }
        else{
            ret->hvalue = ret->root->hash_value();
        }
        return ret;
    }
}

std::string SafraTree::to_string() const{
    using namespace std;
    string ret;
    ret.append("SafraTree #");
    ret.append(INT_TO_STR(this->treeID));
    ret.append(". Name:");
    ret.append(INT_TO_STR(this->name));
    
    ret.append("\n");
    
    string temp;
    boost::to_string(this->used_node_names, temp);
    ret.append("  Used node names: ");
    ret.append( temp );
    ret.append("\n");
    boost::to_string(this->marked_nodes, temp);
    ret.append("  Marked nodes: ");
    ret.append( temp );
    ret.append("\n");
    boost::to_string(this->temp_node_names, temp);
    ret.append("  Temp. reserved names: ");
    ret.append( temp );
    ret.append("\n");

    
    if(this->root != NULL){
        ret.append(this->root->to_string(0));
    } else {
        ret.append("(no nodes)\n");
    }
    
    return ret;
}

std::string SafraTree::to_GASt_string() const{
    if(this->root != NULL){
        return this->root->to_GASt_string(0);
    } else {
        return std::string("(no nodes)\n");
    }
}


SafraTree* SafraTree::get_tree(int i){
    if(!SAVE_TREE_DATA)
        return NULL;
    else
        return SafraTree::trees[i];
}

/*** Implementation of SafraNode ***/

SafraNode::SafraNode(){
    this->marked = false;
}

SafraNode::~SafraNode(){
}

/** Recursive: frees my name and then destructs all of my children. 
 * DO NOT USE for cleaning up an entire tree; just call the tree destructor.
 * Used to pare a branch of nodes and all their descendants from the tree.
 */
void SafraNode::kill_node(){
    this->tree->free_node_name(this->name);
    for(int i = 0; i < this->children.size(); i++)
        this->children[i]->kill_node();
}

bool SafraNode::operator==(const SafraNode& other) const {
    if(this->name != other.name ) return false;
    if( !(this->states == other.states) ) return false;
    if(this->children.size() != other.children.size()){
        return false;
    }
    for(int i = 0; i < this->children.size(); i++){
        if( !( *(this->children[i]) == *(other.children[i])) )
            return false;
    }
    return true;
}

void SafraNode::accumulate_subtree_names(state_set_t& subtree_names){
    subtree_names.set(this->name - 1);
    std::vector<SafraNode*>::iterator iter = this->children.begin();
    while(iter != this->children.end() ){    
        (*iter)->accumulate_subtree_names(subtree_names);
        iter++;
    }
}

/**
 * Updates a Safra tree using one depth-first search.
 * Copy nodes (erasing all marks), create children, and perform transitions of labels. 
 * Then suppress states, suppress node labelings, and mark appropriate nodes for
 * Safra's construction. 
 * @return the new root node of the tree, or NULL if the tree (or subtree) is empty.
 */
SafraNode* SafraNode::clone_spawn_and_transition(SafraTree* new_tree, const NBW& input, int character, bool root, state_set_t& kill_set, const state_set_t& nbw_final_states){    
    // clone this node, transitioning the labels
    void* new_node_loc = (void*)(new_tree->node_storage + (this->name-1));
    
    SafraNode* ret = new(new_node_loc) SafraNode();
    ret->marked = false;
    ret->ID = SafraNode::next_id++;
    ret->name = this->name;
    ret->states = this->states;
    ret->tree = new_tree;
    if(root)
        new_tree->root = ret;

    // if(TRANSITION_FIRST) // Screw this; TRANSITION_FIRST is now mandatory.
    input.transition(ret->states, character);
        
    /** Perform the "eliminate states that my left siblings have, and kill me
     * if I'm empty" steps on the new root node of the subtree.
     */
    if(ret->states.is_subset_of(kill_set)){  // Kill this subtree.
        if(root)
            new_tree->root = NULL;
        else {
            /* This is downright horrible. If this subtree has children, we need
             * to identify all of their names (so that we can mark the names as
             * unused).
             */
            state_set_t names_to_free;
            names_to_free.resize(2*input.size);
            this->accumulate_subtree_names(names_to_free);
            new_tree->temp_node_names |= names_to_free;
            //ret->kill_node();
        }
        return NULL;
    }        
    
    /**
     * OK, this is a tricky bit. We want to create a child before recursing, so
     * that the nodes are named properly. However, we want to add it to the list
     * of children *after* recursing, so that the order of children is preserved
     * correctly in the new, cloned tree. 
     * To solve this problem, we will reserve a node name before recursing, but
     * we will only create the child after recursing, and only if it is not 
     * killed in the process.
     * If the child is killed (all of its states are present in left siblings),
     * we will free the name when exiting the recursion (if this is the root
     * node).
     */
    int new_child_name = new_tree->name_node();

    ret->states -= kill_set;

    // Recursion!
    std::vector<SafraNode*>::iterator iter = this->children.begin();
    
    while(iter != this->children.end() ){
        SafraNode* cloned_child = (*iter)->clone_spawn_and_transition(new_tree,input,character,false, kill_set, nbw_final_states);
        if(cloned_child != NULL)
            ret->children.push_back(cloned_child);
        iter++;
    }    
    
    
    /* Now, we find out if creating a new child is appropriate.
     * Note that the child will survive only if:
     *  - its initial states (states of this node which are final) are nonempty
     *  - states are still nonempty after removing the states in left siblings (kill set)
     *  - this node is not going to be marked (which would kill the children)
     */
    
    state_set_t new_child_states = ret->states;
    new_child_states &= nbw_final_states;
    new_child_states -= kill_set;
    
    kill_set |= new_child_states;
    
    /* Check to see if this node needs marking. Because of recursion, kill_set
     * currently includes all of our descendant's states (including the 
     * hypothetical new child). If our states are a subset of that set, we
     * should be marked (because of the line "ret->states -= kill set" above,
     * we know that ret has no states in the kill set but not in its children.
     * Therefore, if ret->states is a subset of kill_set, then ret->states is
     * equal to the union of the labels of its children.).
     */
    if(ret->states.is_subset_of(kill_set)){
        /* The node does need marking. Kill all its children and mark it.
         */
        ret->marked = true;
        new_tree->mark_node(ret->name);
        int i = 0, s = ret->children.size();
        for(; i < s; i++)
            ret->children[i]->kill_node();
        ret->children.clear();
        
        /* In this case, we're not going to create the new child. Signal the
         * tree to hold that name in reserve until we finish updating, but then
         * to free the name when the update is finished.
         */
        new_tree->temp_name_node(new_child_name);
        
    } else if(new_child_states.any()){ // Check to see if we should create the new child
        // Success! A new child will be created.        
        SafraNode* new_child;
        void* new_child_loc = (void*)(new_tree->node_storage + (new_child_name - 1));        
        new_child = new(new_child_loc) SafraNode();
        new_child->ID = SafraNode::next_id++;
        new_child->name = new_child_name;
        new_child->marked = MARK_NEW_CHILDREN;
        if(MARK_NEW_CHILDREN)
            new_tree->mark_node(new_child->name);
        new_child->states = new_child_states;
        new_child->tree = new_tree;
        ret->children.push_back(new_child);
    } else {
        /* The new child node would have no states, and is immediately
         * deleted in the last step of Safra's construction. Don't create the
         * node, but signal the tree to keep its name reserved until this 
         * transition is complete, so that any future nodes are named properly.
         */
        new_tree->temp_name_node(new_child_name);
    }    
    
    kill_set |= ret->states;     
    
    
    /* If ret is the root node (root = true), the transition of the tree is
     * now complete. Free any node names that have been temporarily reserved.
     */
    if(root){
        new_tree->free_temp_names();
    }
    
    return ret;    
}

std::string SafraNode::to_string(int indent_level) const{
    std::string s, temp;
    s.append("     ");
    for(int i = 1; i < indent_level; i++)
        s.append("     ");
    if (indent_level > 0)
        s.append(" +-> ");    
    s.append("[");
    s.append(INT_TO_STR(this->name));
    s.append("|");  
    
    boost::to_string(this->states, temp);    
    s.append(temp);
    
    s.append("]");
    if(this->marked)
        s.append("!");
    s.append("\n");
    for(int i = 0; i < this->children.size(); i++)
        s.append(this->children[i]->to_string(indent_level+1));
    return s;
}

/** A slightly more verbose form of the to_string method which 
 * conforms more exactly to the GASt syntax.
 */
std::string SafraNode::to_GASt_string(int indent_level) const{
    std::string s;
    s.append("     ");
    for(int i = 1; i < indent_level; i++)
        s.append("     ");
    if (indent_level > 0)
        s.append(" +-> ");    
    s.append("[");
    s.append(INT_TO_STR(this->name));
    s.append("|");  
    
    bool first = true;
    for(int i = 0; i < this->states.size(); i++){
        if(this->states[i]){
            if(!first){
                s.append(",");
            }else{
                first = false;
            }
            s.append(INT_TO_STR(i));
        }
    }
    
    s.append("]");
    if(this->marked)
        s.append("!");
    s.append("\n");
    for(int i = 0; i < this->children.size(); i++)
        s.append(this->children[i]->to_GASt_string(indent_level+1));
    return s;
}

/** Generate a hash value for the node.
 */
std::size_t SafraNode::hash_value() const{
    std::size_t seed = 0; std::string s;
    boost::hash_combine(seed, name);
    int max = this->states.size() - 3;
    std::vector<unsigned long> bit_bucket(this->states.num_blocks());
    std::insert_iterator<std::vector<unsigned long> > i_itr(bit_bucket, bit_bucket.begin());
    boost::to_block_range(this->states, i_itr);
    boost::hash_range(seed, bit_bucket.begin(), bit_bucket.end());
    for(int i = 0; i < this->children.size(); i++)
        boost::hash_combine(seed, this->children[i]->hash_value());
    return seed;
}

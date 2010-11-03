/** DRW.hpp: specifies a deterministic Rabin automaton which
 *  recognizes languages of (one-way) infinite words.
 *
 *  (c) Joe Gershenson, 2009
 */

#pragma once
#ifndef DRW_H
#define DRW_H

#include <map>
#include <string>
#include <vector>

#include "utils.hpp" // for special types
#include "SafraTree.hpp" // to see safra trees for to_GASt_string()

class DRW{
    private:
    
        /* The number of strongly connected components in the automaton.
         * Like all of the Boost-graph-dependent fields, it is not reliable
         * until generated with build_boost_components.
         * If build_boost_components has not yet been called, this field will
         * have the value -1.
         */
        int num_sccs;
    
        /* A Boost Graph Library representation of the transition
         * matrix. The current implementation IGNORES the characters
         * which label each transition. Thus, while it is useful for
         * answering "is the automaton empty", it does not generate
         * counterexamples.
         */
        BoostGraph* t_matrix;
        
        /* The strongly connected components of the graph. Again, this
         * IGNORES the characters labeling each transition.
         */
        std::vector<int>* sccs;
            
        /*
         * Build a version of the transition matrix using the 
         * Boost Graph Library (http://www.boost.org). It is stored in
         * the private field t_matrix and used for graph algorithms such
         * as finding strongly connected components.
         * The current implementation also calculates the strongly connected
         * components.
         */
        void build_boost_components();
        
        
        /*
         * Used to construct the states of the complemented Rabin automaton.
         */
        class CompState{
            int rabin_state;
            bool in_initial_part; // the bool value is true if this state is part of the initial transition system
            state_set_t s1; // finite rabin pairs hit
            state_set_t s2; // infinite rabin pairs hit
            
            int buchi_index; // can be used to cache the index in the buchi automaton
            
            CompState(int state_set_size);
            bool operator==(const CompState& other);
            
            /* Methods: */
            int get_index(std::vector<CompState*>& seen);
            
            friend class DRW; // DRW has no special access privileges otherwise
        };
    
    public:
    
        /**
         * The number of states in the automaton.
         */
        int size;
        
        /**
         * The size of the automaton's alphabet.
         */
        int alphabet_size;

        /**
         * The initial state of the automaton.
         */
        int initial_state;
    
        /**
         * Not yet integrated into all functions, but useful for generating
         * pretty graphs.
         */
        std::string alphabet;
        
        /**
         * This semantics table details what each character of the alphabet *actually* stands for.
         * If an entry is present, it overrides the alphabet when writing graphs.
         * Used the same way as it is in the Buchi files (NBW.cpp / NBW.hpp).
         */
        std::vector<std::string> char_labels;

        std::vector<RabinPair*> pairs;

        /**
         * The transition matrix.
         * State X transitions on character C to the state stored at 
         * transition_matrix[X-1][C-1].
         */
        std::vector<int*> transition_matrix;

        /*
         * Read a Rabin automaton from a text file.
         * If you call to_string() on a Rabin automaton, write it to "foo.txt",
         * and then call parse("foo.txt"), the result will be equivalent to
         * the original automaton.
         */
        static DRW* parse(char* filename);        

        
        /* Return the transition from @param state on @param character.
         * May be useful for data encapuslation in the future (if for example
         * the transition matrix is stored in a different form).
         */
        int transition(int state, int character) const;
        
        /** Generate a printable version of this automaton.
         */
        std::string to_string();
        
        /** Generate a version of this automaton in the GASt output format.
         * This format shows the Safra trees corresponding to each state,
         * so it may be useful for debugging or understanding the 
         * constructions involved. However, it requires maintaining all the
         * trees in memory, so you must set SAVE_TREE_DATA to true in 
         * @file SafraTree.hpp.
         */
        std::string to_GASt_string();
        
        /*
         * Generate a version of the transition graph suitable for rendering
         * using GraphViz (http://www.graphviz.org).
         */
        std::string to_digraph();
        
        /*
         * Determine if the language of the automaton is empty.
         * Assumes that every state is reachable since Safra's construction only builds
         * the reachable part of the Rabin automaton.
         * Uses the following algorithm to check emptiness -- suppose that the Rabin 
         * pairs are given as (FIN, INF).
         *    For each Rabin acceptance pair, see if that pair can be satisfied:
         *      A state in INF with a path to itself which does not pass through
         *      any states in FIN is required to satisfy the pair.
         *    If the pair can be satisfied, the automaton is not empty.
         * Uses BGL implementation of Tarjan's algorithm to identify strongly connected
         * components.
         */
        bool is_empty();
        
        /*
         * Determine if the language of the automaton is universal.
         * Assumes that every state is reachable since Safra's construction only builds
         * the reachable part of the Rabin automaton.
         * Uses the following algorithm to check universality -- suppose that the Rabin 
         * pairs are given as (FIN, INF). We actually look for not-universal: see if 
         * there is a way to fail (not satisfy) any pair. Thus our algorithm is this:
         *    For each Rabin acceptance pair, see if that pair can be failed:
         *      A state in FIN with a path to itself which does not pass through
         *      any states in INF makes it possible to fail the pair.
         *    If no pair can be failed, the automaton is universal.
         * Uses BGL implementation of Tarjan's algorithm to identify strongly connected
         * components.
         */
        bool is_universal();
                
        /* 
         * Test function: print all strongly connected components of the graph.
         */
        void print_components();
        
        /* 
         * Generate and return a BŸchi automaton which accepts the complement 
         * of the language accepted by this automaton.
         */
        NBW* complement();
        
        DRW();
        ~DRW();
};

class RabinPair{
    public:
        state_set_t infinite;
        state_set_t finite;
        
        RabinPair(int size);
        ~RabinPair();
};

#endif

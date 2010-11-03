/** @file NBW.hpp 
 *  Specifies a nondeterministic Buchi automaton which
 *  recognizes languages of (one-way) infinite words.
 *
 *  @author Joe Gershenson
 */

#pragma once
#ifndef NBW_H
#define NBW_H

#include <map>
#include <string>

#include "utils.hpp"
#include "buchi_gen.hpp"
#include "DRW.hpp"
#include "SafraTree.hpp"
#include "logic.hpp"

class NBW{
    private:
        /********************************* pointer fields *********************************/        
        /**
         * The raw transition matrix.
         * 1D array representing 2D array. Transitions are stored in
         *     (state_from-1)*alphabet_size + (char_on)-1
         */
        state_set_t* transition_matrix;
        
        /**
         * Array of state sets, used for caching transition values.
         * Transitions from X on A are stored at
         *  transition_cache[(X->to_ulong()) * alphabet_size + (A-1)]
         */
        state_set_t* transition_cache;


      /****************************** struct / class fields ******************/
        /** 
         * Designates which states have been projected.
         *
         */
        state_set_t projected_tracks;
        
        /**
         * The initial states of the automaton.
         */
        state_set_t initial;
        /**
         * The final (or accepting) states of the automaton.
         */
        state_set_t final;
        
    // Private helper functions
    static NBW* build_helper(const Conjunction& f, Boundary conditions);
    static NBW* build_conjunction_automaton(const Conjunction& f, Boundary conditions);

    public:    
      /********************************* int fields *********************************/

        /** Whether or not to use a cache for the automaton.
         */
        bool use_cache;
        
        /** Whether the automaton has been trimmed since the last time it was
         * modified (if so it contains only accessible + coaccessible states).
         */
        bool trimmed;
        
        int size;
        int alphabet_size;
        int num_transitions;
        
      /******************************** struct/class fields *********************************/
        /**
         * Not yet integrated into all functions, but useful for generating
         * pretty graphs.
         */
        std::string alphabet;

        /**
         * This semantics table details what each character of the alphabet *actually* stands for.
         * If an entry is present, it overrides the alphabet when writing graphs.
         */
        std::vector<std::string> char_labels;

        /**
         * Labels for the states of the automaton. Not used internally (except for graph printouts),
         * currently just used for debugging purposes. May later signify alphabet structure.
         */
        std::vector<std::string> state_labels;


      /********************************* methods *********************************/
        /** Assists in parsing GASt files.
         * Clients should call @parse instead.
         */
        static NBW* parse_from_GASt(std::istream &input, std::string &buffer);
           
        static NBW* parse(char* filename);
        
        static NBW* build_random_automaton( int states, 
            int alphabet_size, 
            double transition_density, 
            double final_state_density );
        
        /** Build an automaton to recognize the specified DNF logical formula,
         * for a cellular automaton with the specified type of boundary 
         * conditions.
         */
        static NBW* build_automaton(std::vector<Conjunction*> formula, Boundary conditions);
       
        static NBW* disjoint_sum(NBW* one, NBW* two);

        static NBW* product(NBW* one, NBW* two);
        
        static NBW* intersection(NBW* one, NBW* two);
       
        /**
         * Transitions the set of states on the given character. If this
         * automaton is in the superposition of states provided, and it reads
         * the character indicated, it transitions to some set of states S. When
         * this function returns, states has the value S.
         * Naturally, this trips an assertion if states not is the size of 
         * this automaton's state set.
         * Now that caching is either done during automaton construction, or 
         * never, this function no longer modifies the NBW object.
         */
        void transition(state_set_t& states, int character) const;
        
        /**
         * Produce a string representation of the automaton. If saved to a file
         * the result may be read by @function parse().
         */
        std::string to_string() const;
        
        /**
         * Produce a string representation of the automaton suitable for 
         * viewing with dot (see http://graphviz.org).
         */        
        std::string to_digraph() const;


        /* Get a copy of the initial states of the automaton.
         * Allowing access to only read-only copies helps for encapsulation.
         */
        state_set_t get_initial_states() const;
        
        /* Get a copy of the final states of the automaton.
         * Allowing access to only read-only copies helps for encapsulation.
         */
        state_set_t get_final_states() const;
        
        /**
         * Returns a bitvector of the accessible states in the automaton.
         */
        state_set_t accessible_states() const;
        
        /**
         * Returns a bitvector of the accessible states in the automaton. A 
         * coaccessible state is one with a path to infinitely many final states.
         */
        state_set_t coaccessible_states() const;


        /** Returns true IFF the language of the automaton is empty. */        
        bool is_empty();
        
        /*
         * Return a deterministic Rabin automaton accepting the same language;
         * uses Safra's construction.
         */
        DRW* determinize() const;        

        /* 
         * Generate and return a Büchi automaton which accepts the complement 
         * of the language accepted by this automaton. Uses determinization.
         * Not const because the automaton is trimmed first.
         */
        NBW* get_complement();
        
        /* 
         * Project a character. ("Erase" the track, or make the given track 
         * irrelevant to the behavior of the automaton.)
         */
        void project(int track_index);

        /* Remove any states which are not (accessible and coacdessible) and
         * condense the automaton accordingly. Has no effect if the size of the
         * automaton is 1, so that empty machines can be represented by a single state.
         * Returns the number of states saved by this process (old size - new size).
         */
        int trim();

        /** Does no setting of fields or allocating of extra memory. For functions 
         * which are going to set all fields manually.
         */
        NBW();
        
        /** Sets all fields directly from the parameters provided. For external
         * classes that have already computed the automaton and thus 
         * "know exactly what they want", or for internal functions which are
         * already building an automaton from an adjacency list.
         */
        NBW(int size, 
            int alphabet_size, 
            std::vector<boost::tuple<int, int,int> > adjacency_list,
            state_set_t initial,
            state_set_t final,
            std::vector<std::string> char_labels = std::vector<std::string>(),
            std::vector<std::string> state_labels = std::vector<std::string>());
            
        ~NBW();
};

#endif

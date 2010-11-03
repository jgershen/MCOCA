#pragma once
#ifndef BUCHI_GEN_H
#define BUCHI_GEN_H

#include <vector>
#include <boost/dynamic_bitset.hpp>

#include "logic.hpp"
#include "utils.hpp"
#include "NBW.hpp"

class BuchiState
{
    private:
        static Conjunction current_formula;
        static int next_state_id;
        static int formula_tracks;
    
        /* Used to search through the states seen so far 
         * (inside @function get_state(), for example)
         * while ignoring the sink state and the initial state, which don't
         * match the "real" states that can be transitioned into.
         */
        static const int FIRST_REAL_STATE = 2;

    public:
        bool accept;
        bool sink;
        
        int state_id;        
        
        slice old_slice;
        slice current_slice;
        
        bitvector neg_lits_sat;
    
        static std::vector<BuchiState*> state_list;    
    
        BuchiState();
        ~BuchiState();
        
        BuchiState* get_target(slice next_slice);
        
        static void cleanup();
        static void initialize(const Conjunction& f);
        
        static BuchiState* get_state(const slice& x, const slice& y,
                                    const bitvector& neg_lits_sat);
        static BuchiState* get_sink();
        
        friend class NBW;
};

#endif
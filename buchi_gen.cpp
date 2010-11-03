/* @file buchi_gen.cpp
 * @author Joe Gershenson
 *
 * Generates Buchi automata to check properties of elementary cellular automata.
 */

#include <vector>

#include <boost/dynamic_bitset.hpp>
#include <boost/tuple/tuple.hpp>

#include "logic.hpp"
#include "buchi_gen.hpp"

//TODO: remove this dependency (it's for debugging)
#include <iostream>

int initial_state_id;

// initialize static members
std::vector<BuchiState*> BuchiState::state_list;
int BuchiState::next_state_id;
Conjunction BuchiState::current_formula; // needs to be reset
int BuchiState::formula_tracks;

/*************** Implementation of BuchiState ***********************/

BuchiState::BuchiState(){
}

BuchiState::~BuchiState(){
}

void BuchiState::cleanup(){
    for(int i = 0; i < state_list.size(); i++)
        delete state_list[i];
    state_list.clear();
}

void BuchiState::initialize(const Conjunction& f){
    BuchiState::formula_tracks = SymbolTable::var_count();

    for(int i = 0; i < state_list.size(); i++)
        delete state_list[i];
    state_list.clear();
    next_state_id = 0;
    BuchiState::current_formula = f;
    BuchiState* sink = new BuchiState();
    sink->state_id = next_state_id++; //TODO(jgershen): fix this to use a single initial state...
    sink->accept = false; sink->sink = true;
    sink->neg_lits_sat = bitvector(f.neg_literals.size());
    state_list.push_back(sink);
    
    BuchiState* initial_state = new BuchiState();
    initial_state->state_id = next_state_id++;
    initial_state->accept = false; initial_state->sink = false;
    initial_state->neg_lits_sat = bitvector(f.neg_literals.size());
    state_list.push_back(initial_state);
    initial_state_id = initial_state->state_id;
    
}

BuchiState* BuchiState::get_sink(){
    return state_list[0];
}

BuchiState* BuchiState::get_state(const slice& x, const slice& y,
                                    const bitvector& neg_lits_sat){
    /* look for a previously created state with these attributes; ignore the sink
     * and the initial state
     */
    for(int i = FIRST_REAL_STATE; i < state_list.size(); i++){
        if( state_list[i]->old_slice == x 
            && state_list[i]->current_slice == y 
            && state_list[i]->neg_lits_sat == neg_lits_sat)
            return state_list[i];    
    }
    BuchiState* new_state = new BuchiState();
    new_state->sink = false;
    new_state->state_id = next_state_id++;
    new_state->old_slice = x;
    new_state->current_slice = y;
    new_state->neg_lits_sat = neg_lits_sat;
    new_state->accept = (neg_lits_sat.count() == current_formula.neg_literals.size());
    std::string s1, s2;
    boost::to_string(x,s1); boost::to_string(y,s2);
    state_list.push_back(new_state);
    
    return new_state;
}

BuchiState* BuchiState::get_target(slice next_slice){
    bool pos_lits_ok = current_formula.check(this->old_slice, this->current_slice, next_slice);
    
    if(pos_lits_ok){
        bitvector neg_lits_sat = bitvector(this->neg_lits_sat);
        for(int i = 0; i < current_formula.neg_literals.size(); i++){
            if(! current_formula.neg_literals[i]->check(this->old_slice, this->current_slice, next_slice))
                neg_lits_sat.set(i);
        }
        return BuchiState::get_state(this->current_slice, next_slice, neg_lits_sat);
    }else{   
        return BuchiState::get_sink();
    }
}



NBW* NBW::build_automaton(std::vector<Conjunction*> formula, Boundary conditions){
    NBW* ret = build_conjunction_automaton(*formula[0], conditions);
    for(int i = 1; i < formula.size(); i++){
        NBW* tmp = build_conjunction_automaton(*formula[i], conditions);
        NBW* sum = NBW::disjoint_sum(ret,tmp);
        delete ret;
        delete tmp;
        ret = sum;
    }
    return ret;
}


/* Does the parts of building an NBA which don't involve quantifiers. Does not
 * handle negative literals, either.
 */
NBW* NBW::build_helper(const Conjunction& f, Boundary conditions){
    BuchiState::initialize(f);

    /* 
     * Generate all characters in the alphabet and put them in a vector.
     */
    std::vector<slice> characters;
    for(unsigned long c = 0; c < (1 << (BuchiState::formula_tracks)); c++){
        characters.push_back( slice(BuchiState::formula_tracks, c) );
    }
    
    // Adjacency list format is (state, character, state)
    std::vector<boost::tuple<int, int, int> > adjacency_list;
    std::vector<BuchiState*> work_queue; 
    slice zeroes(BuchiState::formula_tracks, 0);
    slice neg_lits_unsat(f.neg_literals.size(), 0);
    
    if(conditions == OMEGA){
        for(int i = 0; i < characters.size(); i++){
            BuchiState* b = BuchiState::get_state(zeroes, characters[i], neg_lits_unsat);
            work_queue.push_back(b);
            adjacency_list.push_back(boost::make_tuple(initial_state_id, i, b->state_id));
        }    
    } else if (conditions == ZETA) {
        for(int i = 0; i < characters.size(); i++){
            for(int j = 0; j < characters.size(); j++){
                BuchiState* b = BuchiState::get_state(characters[i], characters[j], neg_lits_unsat);
                work_queue.push_back(b);
                adjacency_list.push_back(boost::make_tuple(initial_state_id, j, b->state_id));
            }
        }       
    }
    
    while( !work_queue.empty() ){
        std::vector<BuchiState*> new_work_queue;
        
        int max = work_queue.size();

        // calculate Buchi transitions
        for(int i = 0; i < max; i++){
            for(int j = 0; j < characters.size(); j++){
                int x = BuchiState::next_state_id;                    
                BuchiState* new_state = work_queue[i]->get_target(characters[j]);
                adjacency_list.push_back(boost::make_tuple(work_queue[i]->state_id, j, new_state->state_id));
                if(x < BuchiState::next_state_id){ // state is new
                    new_work_queue.push_back(new_state);
                }
            }
        }
        
        work_queue = new_work_queue;
    }
    
    for(int i = 0; i < BuchiState::state_list.size(); i++)
        assert(BuchiState::state_list[i]->state_id == i);

    /* We have finished calculating the transition matrix; build the automaton. */
            
    int size = BuchiState::next_state_id;
    int alphabet_size = characters.size();

    std::vector<std::string> char_labels;
    // store character semantics table
    for(int i = 0; i < characters.size(); i++){
        std::string s; 
        boost::to_string(characters[i], s);
        char_labels.push_back(s);
    }
    
    std::vector<std::string> state_labels;
    // store state semantics table (state labels);
    state_labels.push_back("SINK");
    state_labels.push_back("INITIAL");
    for(int i = BuchiState::FIRST_REAL_STATE; i < size; i++){
        assert(i == BuchiState::state_list[i]->state_id);
        std::string s1; std::string s2;
        boost::to_string(BuchiState::state_list[i]->old_slice, s1);
        boost::to_string(BuchiState::state_list[i]->current_slice, s2);
        s1 += ":";
        s1 += s2;        
        state_labels.push_back(s1);    
    } 
    
    // set initial and final states
    state_set_t initial(size);
    initial.set(initial_state_id);
    state_set_t final(size);
    for(int i = 0; i < size; i++){
        if(BuchiState::state_list[i]->accept)
            final.set(i);
    }
    
    NBW* ret = new NBW(size, 
                   alphabet_size, 
                   adjacency_list,
                   initial,
                   final,
                   char_labels,
                   state_labels);

    // Free memory & data structures & things
    BuchiState::cleanup();
    
    return ret;
}



NBW* NBW::build_conjunction_automaton(const Conjunction& f, Boundary conditions){        
    NBW* ret = build_helper(f, conditions);
    
    bool current_formula_needs_negated = false; // for double negatives
    
    // Handle quantifiers here
    for(int i = f.quantifiers.size() - 1; i >= 0; i--){
        if(f.quantifiers[i].universal){ // universal quantifiers w/negs
            if(f.quantifiers[i].negated){
                if(current_formula_needs_negated){
                    // ~Ax ~foo --> Ex foo
                    current_formula_needs_negated = false;
                    ret->project(f.quantifiers[i].variable_index);                
                } else {
                    // ~Ax foo -> Ex ~foo
                    current_formula_needs_negated = false;                
                    NBW* not_foo = ret->get_complement();
                    delete ret;
                    ret = not_foo;
                    ret->project(f.quantifiers[i].variable_index);
                }
            } else { // universal quantifiers w/o negation to immediate left
                if(current_formula_needs_negated){
                    // Ax ~foo -> ~Ex foo
                    ret->project(f.quantifiers[i].variable_index);
                    current_formula_needs_negated = true;                    
                    // still needs negated
                } else {
                    // Ax foo --> ~Ex ~foo
                    current_formula_needs_negated = true;
                    NBW* not_foo = ret->get_complement();
                    delete ret;
                    ret = not_foo;
                    ret->project(f.quantifiers[i].variable_index);                    
                }
            }
        } else { // existential quantifiers w/negs
            if(f.quantifiers[i].negated){
                if(current_formula_needs_negated){
                    // ~Ex ~foo
                    NBW* not_foo = ret->get_complement();
                    delete ret;
                    ret = not_foo;
                    ret->project(f.quantifiers[i].variable_index);                  
                    current_formula_needs_negated = true;
                } else {
                    // ~Ex foo
                    ret->project(f.quantifiers[i].variable_index);
                    current_formula_needs_negated = true;
                }
            } else { // existential quantifiers w/o negs
                if(current_formula_needs_negated){
                    // Ex ~foo
                    NBW* not_foo = ret->get_complement();
                    delete ret;
                    ret = not_foo;
                    ret->project(f.quantifiers[i].variable_index);                  
                    current_formula_needs_negated = false;
                } else {
                    // Ex foo
                    ret->project(f.quantifiers[i].variable_index);
                    current_formula_needs_negated = false;
                }
            }
        }
        
        ret->trim(); // trim between each cycle
    } // move on to next quantifier
    
    // negation propogated to topmost level
    if(current_formula_needs_negated){
        NBW* not_foo = ret->get_complement();
        delete ret;
        ret = not_foo;
        ret->trim();
    }
    
    return ret;
}

/** @file NBW.hpp 
 *  Implementation of a nondeterministic Buchi automaton which
 *  recognizes languages of (one-way) infinite words.
 *  
 *  For specification, see @file NBW.hpp.
 *  @author Joe Gershenson
 */

#include "NBW.hpp"
#include "utils.hpp"

/* Include string IO + stream utilities for tokenizing input files */
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>

/* OpenMP library for parallel operations */
#include <omp.h>

#include <iomanip> // TODO: remove; used for debugging ouput only

#include <boost/unordered_set.hpp>
#include <boost/graph/strong_components.hpp>
#include <boost/graph/reverse_graph.hpp>

#include <boost/graph/adjacency_iterator.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/graph_utility.hpp>

#include "SafraTree.hpp"

/**
 * Hash set type used for hashing SafraTrees.
 */
typedef boost::unordered_set<SafraTree*, stp_hash_t, stp_eq_t> stree_set_t;

NBW* NBW::parse_from_GASt(std::istream &input, std::string &buffer){
    using namespace std;
    NBW* ret = new NBW();
    ret->size = atoi(buffer.c_str()); //read size (already stored in buffer)
    ret->use_cache = (NBW_USE_CACHE && ret->size <= NBW_MAX_CACHED_SIZE);
    
    get_next_line(input, buffer);
    ret->alphabet = buffer;
    
    ret->alphabet_size = buffer.size();
    
    ret->projected_tracks.resize(ret->alphabet_size);    

    // A GASt format automaton always uses the first state as the initial state
    ret->initial.resize(ret->size);
    ret->initial.set(0);
    
    /* Read final states -- in the examples I have, there is only one final
     * state for a given automaton. It's unclear how formatting would change
     * for more than one final state, so I only attempt to read one final state.
     */
    get_next_line(input, buffer);
    ret->final.resize(ret->size);
    ret->final.set(atoi(buffer.c_str()));
    
    // prepare to read transition matrix
    ret->num_transitions = 0;
    
    // allocate memory for transition matrix
    ret->transition_matrix = new state_set_t[ (ret->alphabet_size) * (ret->size)];
    for(int i = 0; i < (ret->alphabet_size)*(ret->size); i++){
        ret->transition_matrix[i].resize(ret->size);
    }

    get_next_line(input, buffer);
    // process transitions
    do{
        string part;
        istringstream iss(buffer);
        getline(iss, part, ' ');
        int state_from = atoi(part.c_str());
        getline(iss, part, ' '); 
        int char_on = 0;
        assert(part.length() == 1);
        while(part[0] != ret->alphabet[char_on])
            char_on++;
        getline(iss, part, ' ');
        int state_to = atoi(part.c_str());
        int index = (state_from)*(ret->alphabet_size)  + (char_on);
        ret->transition_matrix[index].set(state_to);
        ret->num_transitions++;
        get_next_line(input, buffer);
    } while ( !input.eof() );

    
    if(ret->use_cache){
        int cache_size = (ret->alphabet_size) * (1<<(ret->size));
    
        // allocate memory for transition cache
        ret->transition_cache = new state_set_t[cache_size];

        /* Build the entire cache.
         */
        for(unsigned long states = 0; states < (1 << (ret->size)); states++){
            state_set_t states_from(ret->size, states);
            for(int c = 0; c < ret->alphabet_size; c++){   
                int cache_index = (states*(ret->alphabet_size)) + c;
                ret->transition_cache[cache_index].resize( ret->size );
                for(int i = 0; i < ret->size; i++){
                    if(states_from[i])
                        ret->transition_cache[cache_index] |= ret->transition_matrix[ i * ret->alphabet_size + c];
                }
            }
        }
        
    }

    return ret;
}

static int compare_stateset_pointers(state_set_t* one, state_set_t* two){
    return (*one < *two);
}

NBW* NBW::parse(char* filename) {
    using namespace std;
    string s;

    ifstream inf (filename);
    if (inf.good()) {
        get_next_line(inf, s);


        // check to see if file is in GASt format
        if( (s.compare("BUECHI") != 0) && (s.compare("BUCHI") != 0) ){
            return parse_from_GASt(inf, s);
        }
        
        // allocate memory
        NBW* ret = new NBW();

        // determine automaton size (number of states)
        get_next_line(inf, s);
        ret->size = atoi(s.c_str());
        ret->use_cache = (NBW_USE_CACHE && ret->size <= NBW_MAX_CACHED_SIZE);        

        // determine alphabet size
        get_next_line(inf, s);
        ret->alphabet_size = atoi(s.c_str());

        // get number of transitions (old stuff, but part of the input file...)
        get_next_line(inf, s);	
        ret->num_transitions = atoi(s.c_str());
        
        // allocate memory for transition matrix
        ret->transition_matrix = new state_set_t[ (ret->alphabet_size) * (ret->size)];
        for(int i = 0; i < (ret->alphabet_size)*(ret->size); i++){
            ret->transition_matrix[i].resize(ret->size);
        }        

        // process transitions
        for(int i = 0; i < ret->num_transitions; i++){
            get_next_line(inf, s);
            string part;
            istringstream iss(s);
            getline(iss, part, ' ');
            int state_from = atoi(part.c_str());
            getline(iss, part, ' ');
            getline(iss, part, ' ');
            int char_on = atoi(part.c_str());
            getline(iss, part, ' ');
            getline(iss, part, ' ');
            int state_to = atoi(part.c_str());
            int index = (state_from-1)*(ret->alphabet_size)  + (char_on - 1);
            ret->transition_matrix[index].set(state_to - 1);        
        }
        
        if(ret->use_cache){
            int cache_size = (ret->alphabet_size) * (1<<(ret->size));
        
            // allocate memory for transition cache
            ret->transition_cache = new state_set_t[cache_size];
    
            /* Build the entire cache.
             */
            for(unsigned long states = 0; states < (1 << (ret->size)); states++){
                state_set_t states_from(ret->size, states);
                for(int c = 0; c < ret->alphabet_size; c++){   
                    int cache_index = (states*(ret->alphabet_size)) + c;
                    ret->transition_cache[cache_index].resize( ret->size );
                    for(int i = 0; i < ret->size; i++){
                        if(states_from[i])
                            ret->transition_cache[cache_index] |= ret->transition_matrix[ i * ret->alphabet_size + c];
                    }
                }
            }
            
        }


        ret->projected_tracks.resize(ret->alphabet_size);

        // allocate memory for initial state set
        ret->initial.resize(ret->size);
        // assign initial state(s)
        get_next_line(inf, s);        
        string part;
        istringstream iss(s);
        getline(iss, part, ' ');
        do{
            int initial = atoi(part.c_str());
            ret->initial.set(initial - 1);
            getline(iss, part, ' ');
        } while (!iss.eof() && part.size() > 0);
                
        // allocate memory for final state set
        ret->final.resize(ret->size);
        // assign final state(s)
        get_next_line(inf, s);        
        istringstream iss2(s);
        getline(iss2, part, ' ');       
        while(!iss2.eof() && part.size() > 0){
            int final = atoi(part.c_str());
            ret->final.set(final - 1);
            getline(iss2, part, ' ');       
        }
        
        return ret;
    }
    // code reaching this point implies IO error.
    std::cout << "I/O error with " << filename << "; aborting.\n";
    return 0;
}

std::string NBW::to_digraph() const{
    using namespace std;
    bool using_state_labels = state_labels.size() > 0;
    bool using_char_labels = char_labels.size() > 0;
    
    string s;
    s.reserve(300);    
    s.append("digraph buchi_automaton { \n");

    s.append("node [shape=circle];\n");
    
    // mark final states as final, and output state labels for all states if used
    for(int i = 0; i < this->size; i++){
        if(final[i]){
            s.append(INT_TO_STR(i+1));
            s.append(" [peripheries=2");
            if(using_state_labels){
                s.append(",label=\"");
                s.append(state_labels[i]);
                s.append("\"");
            }
            s.append("];\n");
        } else {
            s.append(INT_TO_STR(i+1));
            if(using_state_labels){
                s.append(" [label=\"");
                s.append(state_labels[i]);
                s.append("\"]");
            }
            s.append(";\n");
        }
    }
    
    // mark initial states as initial, using in-arrows from invisible states
    for(int i = 0; i < this->size; i++){
        if(initial[i]){
            s.append("I");
            s.append(INT_TO_STR(i+1));
            s.append(" [style=invis];\n");            
            s.append("I");
            s.append(INT_TO_STR(i+1));
            s.append(" -> ");
            s.append(INT_TO_STR(i+1));
            s.append(";\n");
        }            
    }
 
    for(int s1 = 0; s1 < this->size; s1++){
        for(int s2 = 0; s2 < this->size; s2++){
            string label;
            for(int c = 0; c < this->alphabet_size; c++){
                if(this->transition_matrix[s1*(this->alphabet_size)+c][s2]){
                    if(label.size() > 0)
                        label.append(",");
                    if(using_char_labels)
                        label.append( char_labels[c] );
                    else
                        label.append( INT_TO_STR(c+1) );
                }
            }            
            if(label.size() > 0){
                s.append( INT_TO_STR(s1+1) );
                s.append( " -> " );
                s.append( INT_TO_STR(s2+1) );
                s.append(" [label=\"");
                s.append(label);
                s.append("\"]");
                s.append( ";\n" );
            }
    	}
    }

    s.append("}\n");
    return s;
}

std::string NBW::to_string() const{
    using namespace std;
    string s;
    s.reserve(300);    
    s.append("#----- omega-automaton (NBW) ----- \n");
    s.append("BUECHI\n");
    s.append("# Number of states: \n");
    s.append(INT_TO_STR(this->size));

    //write state semantics if present
    for(int i = 0; i < this->state_labels.size(); i++){
        s.append("\n# ");
        s.append(INT_TO_STR(i+1));
        s.append(":");
        s.append(this->state_labels[i]);
    }

    s.append("\n# Size of alphabet: \n");
    s.append(INT_TO_STR(this->alphabet_size));
    
    //write character labels if present
    for(int i = 0; i < this->char_labels.size(); i++){
        s.append("\n# ");
        s.append(INT_TO_STR(i+1));
        s.append(":");
        s.append(this->char_labels[i]);
    }
    
    s.append("\n# Number of transitions: \n");
    s.append(INT_TO_STR(this->num_transitions));
    s.append("\n# List of transitions: \n");
    for(int i = 0; i < this->size; i++){
        for(int j = 0; j < this->alphabet_size; j++){
            state_set_t targets = this->transition_matrix[i*(this->alphabet_size)  + j];
            for(int k = 0; k < this->size; k++){
                if(targets[k]){
                    s.append( INT_TO_STR(i+1) );
                    s.append( " > " );
                    s.append( INT_TO_STR(j+1) );
                    s.append( " > " );
                    s.append( INT_TO_STR(k+1) );
                    s.append( "\n" );
                }
    	    }
    	}
    }
    s.append("# Initial state(s)\n");
    for(int i = 0; i < this->size; i++){
        if(this->initial[i]){
            s.append(INT_TO_STR(i+1));
            s.append("\n");
        }
    }
    s.append("# Final state(s)\n");
    for(int i = 0; i < this->size; i++){
        if(this->final[i]){
            s.append(INT_TO_STR(i+1));
            s.append(" ");
        }
    }
    s.append("\n");
    s.append("# EOF\n");
    return s;
}

state_set_t NBW::get_initial_states() const{
    return state_set_t(this->initial);   
}

/* Get a copy of the final states of the automaton.
 */
state_set_t NBW::get_final_states() const{
    return state_set_t(this->final);   
}

NBW::NBW(){
    this->trimmed = false;
}

/** Sets all fields directly from the parameters provided. For external
 * classes that have already computed the automaton and thus 
 * "know exactly what they want", or for internal functions which are
 * already building an automaton from an adjacency list.
 */
NBW::NBW(int size, 
        int alphabet_size, 
        std::vector<boost::tuple<int, int,int> > adjacency_list,
        state_set_t initial,
        state_set_t final,
        std::vector<std::string> char_labels,
        std::vector<std::string> state_labels){
    
    this->size = size;
    this->alphabet_size = alphabet_size;
    this->alphabet = default_alphabet;
    this->trimmed = false;
    this->initial = initial;
    this->final = final;
    this->num_transitions = adjacency_list.size();
    if(char_labels.size() > 0) 
        this->char_labels = std::vector<std::string>(char_labels);
    if(state_labels.size() > 0) 
        this->state_labels = std::vector<std::string>(state_labels);
    
    // allocate memory for transition matrix
    this->transition_matrix = new state_set_t[ (this->alphabet_size) * (this->size)];
    for(int i = 0; i < (this->alphabet_size)*(this->size); i++){
        this->transition_matrix[i].resize(this->size);
    }
    
    for(int i = 0; i < adjacency_list.size(); i++){
        int s1 = adjacency_list[i].get<0>();
        int c = adjacency_list[i].get<1>();
        int s2 = adjacency_list[i].get<2>();
        
        this->transition_matrix[ (s1 * this->alphabet_size) + c ].set(s2);
    }       
        
    this->use_cache = (NBW_USE_CACHE && this->size <= NBW_MAX_CACHED_SIZE);            
    // Cache the Buchi automaton's transitions, if it is small enough...
    if(this->use_cache){
        int cache_size = (this->alphabet_size) * (1<<(this->size));    
        // Allocate memory for transition cache
        this->transition_cache = new state_set_t[cache_size];
        // Build the cache.
        for(unsigned long states = 0; states < (1 << (this->size)); states++){
            state_set_t states_from(this->size, states);
            for(int c = 0; c < this->alphabet_size; c++){   
                int cache_index = (states*(this->alphabet_size)) + c;
                this->transition_cache[cache_index].resize( this->size );
                for(int i = 0; i < this->size; i++){
                    if(states_from[i])
                        this->transition_cache[cache_index] |= this->transition_matrix[ i * this->alphabet_size + c];
                }
            }
        }        
    }    
    
} // end NBW(...) -- transition list constructor

NBW::~NBW(){
    // delete the transition matrix
    delete [] this->transition_matrix;
    
    // delete the transition cache
    if (this->use_cache)
        delete [] this->transition_cache;
        
    SafraTree::reset();
}

/** Transition the given set of states on the given character. Alters the calling arguments.
 * Internally, attempts to look up the state set in the cache. If it finds a mapping, 
 * uses that mapping to set the value. If it does not, it calculates the value, caches the mapping
 * (creating new instances for 'from' and 'to' to do so) and then sets the value. In this case
 * two constructor calls are made.
 * Accordingly, the destructor for NBW deletes the cache.
 *
 */
void NBW::transition(state_set_t& states_from, int character) const{
    //assert(character > 0); // we're doing the subtract-by-one internally here.
    if(this->use_cache){
        int cache_index = (states_from.to_ulong() * this->alphabet_size) + (character-1);
        states_from = this->transition_cache[cache_index];
    } else {
        state_set_t temp(this->size);
        for(int i = 0; i < this->size; i++){
            if(states_from[i])
                temp |= this->transition_matrix[ i * this->alphabet_size + (character-1)];
        }
        states_from = temp;
    }
}

DRW* NBW::determinize() const{
    using namespace std;
    
    SafraTree::reset();
    vector<SafraTree*> work_queue;
    vector<SafraTree*> tree_list; //actual states of the automaton       
    stree_set_t trees;

    int states_seen = 1;


    // Create Rabin automaton.
    DRW* ret = new DRW();
    ret->alphabet = this->alphabet;
    ret->char_labels = this->char_labels;
    ret->alphabet_size = this->alphabet_size;
    ret->initial_state = 0;

    /* If the Buechi automaton (this) is empty, then the Rabin
     * automaton is also empty and there is no point in determinizing.
     * To keep things standard, though, we will return a 1-state
     * Rabin automaton with no transitions.
     */ 
    if(this->size == 0){
        ret->size == 1;
        // TODO: add transitions and an acceptance pair to fill empty automaton
        return ret;
    }

    ret->transition_matrix.push_back(new int[ret->alphabet_size]);

    //double stime = 0.0;
    //double ptime = 0.0;

    
    // create initial state
    SafraTree* initial_state = SafraTree::build_initial_tree(*this);
    trees.insert(initial_state);
    tree_list.push_back(initial_state);
    
    // add initial work units
    for(int i = 1; i <= this->alphabet_size; i++)
        work_queue.push_back(initial_state);
    
    while( !work_queue.empty() ){
        if( ((work_queue.size() + tree_list.size()) * this->alphabet_size) > 10000){
            int x = tree_list.size();
            int y = work_queue.size();
            float per = (float(x)) / float(x+y) * 100.0;
            printf("done: %d todo: %d (%.4f%%)\n", x, y,per);
        }
        
        // cout << "work_queue size: " << INT_TO_STR( work_queue.size() ) << endl;
        
        vector<SafraTree*> new_work_queue;
        
        int max = work_queue.size();

        //double parallel_start = omp_get_wtime();

        // ----------- TRANSITION -----------
        //#pragma omp parallel for
        for(int i = 0; i < max; i++){
            for(int j = 0; j < this->alphabet_size; j++){
            // perform transition
                work_queue[i]->targets[j] = SafraTree::get_transition(*work_queue[i], *this, j+1);
            }
        }

        //double parallel_end = omp_get_wtime();
        //ptime += (parallel_end - parallel_start); 
        
//        #pragma omp barrier
        //double serial_start = omp_get_wtime();
        
        // --- reduce (name states) ---
        for(int i = 0; i < max; i++){
            for(int j = 0; j < this->alphabet_size; j++){
                SafraTree* result = work_queue[i]->targets[j];
                
                stree_set_t::iterator loc = trees.find(result);            
                
                if(loc == trees.end()){ // not found
                    result->name = states_seen++;
                    tree_list.push_back(result);
                    trees.insert(loc, result);
                    new_work_queue.push_back(result);
                } else if (result->treeID != (*loc)->treeID) { // tree already reached.
                    work_queue[i]->targets[j] = *loc; //canonicalize reference
                }  
            }
        }
                
        //work_queue.clear();
        work_queue = new_work_queue;
        //new_work_queue.clear();
        
        //double serial_end = omp_get_wtime();
        //stime += (serial_end - serial_start); 
    }
    
    //cout << "   Parallel region used " << setiosflags(ios::fixed) << setprecision(4) << ptime << " s." << endl;
    //cout << "   Serial region used " << setiosflags(ios::fixed) << setprecision(4) << stime << " s." << endl;
    
    ret->size = trees.size();


    /* convert the network of trees to a transition matrix */
    int max = tree_list.size();
    for(int i = 0; i < max; i++){
        ret->transition_matrix.push_back(new int[ret->alphabet_size]);
    }
    
    //#pragma omp parallel for
    for(int i = 0; i < max; i++){
        for(int j = 0; j < ret->alphabet_size; j++)
            ret->transition_matrix[i][j] = tree_list[i]->targets[j]->name;
    }
    
    // THERE IS A LOCK HERE, do not comment it out by accident!
    // jesus christ people.
    //#pragma omp parallel for
    for(int i = 0; i < 2 * (this->size); i++){
        RabinPair* pair = new RabinPair(ret->size);
        for(int j = 0; j < ret->size; j++){
            if(tree_list[j]->marked_nodes[i])
                pair->infinite.set(j);
            else if (tree_list[j]->used_node_names[i] == false)
                pair->finite.set(j);
        }
        if(pair->infinite.any()){
            //#pragma omp critical (ret_pair_locks) // OMG LOCK.
                ret->pairs.push_back(pair);
        } else {
            delete pair;
        }
    }


    SafraTree::trees = tree_list;

    /* If SAVE_TREE_DATA is set to true, we save the data from the Safra tree
     * construction for debugging/output purposes.
     * Be sure to free it later by calling SafraTree::reset!
     */
    if(!SAVE_TREE_DATA){
        SafraTree::reset();
    }

    return ret;
} // end DRW* NBW::determinize() const


NBW* NBW::build_random_automaton( int states, int alphabet_size, double transition_density, double final_state_density ){
    NBW* ret = new NBW();
    ret->size = states; //read size (already stored in buffer)
    ret->use_cache = (NBW_USE_CACHE && ret->size <= NBW_MAX_CACHED_SIZE);
    
    ret->alphabet = default_alphabet; // defined in utils.cpp
    
    ret->alphabet_size = alphabet_size;
        
    // We will always designate 0 ("state 1", in a string representation) 
    // as the lone initial state for a randomly generated automaton.
    ret->initial.resize(ret->size);
    ret->initial.set(0);
    
    /* Randomly designate final states. Along with (?) Vardi we set the initial
     * state to be final, and each other state has a probability of "final_state_density"
     * of being final.
     */
    ret->final.resize(ret->size);
    
    ret->final.set(0);
    
    for(int i = 1; i < ret->size; i++){
        double num = static_cast<double>( rand() ) / static_cast<double>( RAND_MAX );
        if(num < final_state_density)
            ret->final.set(i);
    }
        
    /* Prepare to assign transitions.
     */
    ret->num_transitions = 0;
    
    // allocate memory for transition matrix
    ret->transition_matrix = new state_set_t[ (ret->alphabet_size) * (ret->size)];
    for(int i = 0; i < (ret->alphabet_size)*(ret->size); i++){
        ret->transition_matrix[i].resize(ret->size);
    }

    /* Randomly activate transitions. Note that in this model any transition 
     * (from state s, on character c, to state s') has a uniform and independent
     * probability of being present or not, so the number of transitions in the
     * automaton will be have a binomial distribution.
     */
    for(int c = 0; c < ret->alphabet_size; c++){
        for(int s = 0; s < ret->size; s++){
            for(int s2 = 0; s2 < ret->size; s2++){
                double num = static_cast<double>( rand() ) / static_cast<double>( RAND_MAX );
                if(num < final_state_density){
                    int index = (s)*(ret->alphabet_size) + (c);
                    ret->transition_matrix[index].set(s2);                
                }
            }
        }
    }

    
    /* Build the cache, if necessary.
     */    
    if(ret->use_cache){
        int cache_size = (ret->alphabet_size) * (1<<(ret->size));
    
        // allocate memory for transition cache
        ret->transition_cache = new state_set_t[cache_size];

        /* Build the entire cache.
         */
        for(unsigned long states = 0; states < (1 << (ret->size)); states++){
            state_set_t states_from(ret->size, states);
            for(int c = 0; c < ret->alphabet_size; c++){   
                int cache_index = (states*(ret->alphabet_size)) + c;
                ret->transition_cache[cache_index].resize( ret->size );
                for(int i = 0; i < ret->size; i++){
                    if(states_from[i])
                        ret->transition_cache[cache_index] |= ret->transition_matrix[ i * ret->alphabet_size + c];
                }
            }
        }
        
    }

    return ret;
} // end NBW* NBW::build_random_automaton( int states, int alphabet_size, double transition_density, double final_state_density )

NBW* NBW::get_complement() {
    this->trim();
    DRW* det = this->determinize();
    NBW* ret = det->complement();
    delete det;
    return ret;
} // end NBW* NBW::complement() const

/**
 * Create the NBW automaton accepting L(one) \union L(two). 
 * It has |one| + |two| states.
 * *Requires that one and two have the same alphabet.*
 */
NBW* NBW::disjoint_sum(NBW* one, NBW* two){
    one->trim();
    two->trim();

    NBW* ret = new NBW();
    ret->size = one->size + two->size;
    ret->use_cache = (NBW_USE_CACHE && ret->size <= NBW_MAX_CACHED_SIZE);
    
    ret->alphabet = one->alphabet;
    ret->alphabet_size = one->alphabet_size;
    
    ret->char_labels = std::vector<std::string>(one->char_labels);
    
    // preserve pretty state labels if both input automata had them
    if(one->state_labels.size() >= one->size 
                            && two->state_labels.size() >= two->size){
        ret->state_labels = std::vector<std::string>();
        for(int i = 0; i < one->size; i++)
            ret->state_labels.push_back("1-" + one->state_labels[i]);
        for(int i = 0; i < two->size; i++)
            ret->state_labels.push_back("2-" + two->state_labels[i]);          
    }
    
    /* We allocate memory in which to store a bitvector of projected tracks --
     * although no tracks in this automaton should be projected yet!
     */
    ret->projected_tracks.resize(ret->alphabet_size);    

    /* The initial states of both one and two are initial states in the disjoint
     * sum automaton, although the initial states of two are relabeled.
     */
    ret->initial.resize(ret->size);
    
    state_set_t tmp = one->initial;
    tmp.resize(ret->size);
    ret->initial = tmp;
    
    for(int i = 0; i < two->size; i++){
        if(two->initial[i])
            ret->initial.set(i + one->size);
    }

    /* We are accepting the union of two languages, so a state accepts if
     * it accepts in either language.
     */
    ret->final.resize(ret->size);
    for(int i = 0; i < one->size; i++){
        if(one->final[i]) ret->final.set(i);
    }
    for(int i = 0; i < two->size; i++){
        if(two->final[i]) ret->final.set(one->size + i);
    }

    /* Set number of transitions */
    ret->num_transitions = one->num_transitions + two->num_transitions;
    
    /* allocate memory for transition matrix */
    ret->transition_matrix = new state_set_t[ (ret->alphabet_size) * (ret->size)];
    for(int i = 0; i < (ret->alphabet_size)*(ret->size); i++){
        ret->transition_matrix[i].resize(ret->size);
    }

    /* import transitions from NBW one */
    for(int s0 = 0; s0 < one->size; s0++){
        for(int c = 0; c < one->alphabet_size; c++){
            state_set_t tmp = one->transition_matrix[s0*one->alphabet_size + c];
            tmp.resize(ret->size);
            ret->transition_matrix[s0*ret->alphabet_size + c] = tmp;
        }
    }

    /* import transitions from NBW two */
    for(int s0 = 0; s0 < two->size; s0++){
        for(int c = 0; c < two->alphabet_size; c++){
            state_set_t tmp = two->transition_matrix[s0*two->alphabet_size + c];
            int retindex = (s0 + one->size)*ret->alphabet_size + c;
            for(int s1 = 0; s1 < two->size; s1++){
                if(tmp[s1])
                    ret->transition_matrix[retindex].set(s1+one->size);
            }
        }
    }

    /* build a cache, if necessary */
    if(ret->use_cache){
        int cache_size = (ret->alphabet_size) * (1<<(ret->size));
        // allocate memory for transition cache
        ret->transition_cache = new state_set_t[cache_size];
        /* Build the entire cache.
         */
        for(unsigned long states = 0; states < (1 << (ret->size)); states++){
            state_set_t states_from(ret->size, states);
            for(int c = 0; c < ret->alphabet_size; c++){   
                int cache_index = (states*(ret->alphabet_size)) + c;
                ret->transition_cache[cache_index].resize( ret->size );
                for(int i = 0; i < ret->size; i++){
                    if(states_from[i])
                        ret->transition_cache[cache_index] |= ret->transition_matrix[ i * ret->alphabet_size + c];
                }
            }
        }
        
    }
    
    return ret;    
} // end NBW* NBW::disjoint_sum(NBW* one, NBW* two)

/**
 * Create the NBW automaton accepting L(one) \product L(two). 
 * It has |one| * |two| states.
 * *Requires that one and two have the same alphabet.*
 */
NBW* NBW::product(NBW* one, NBW* two){
    one->trim();
    two->trim();


    NBW* ret = new NBW();
    ret->size = one->size * two->size;
    ret->use_cache = (NBW_USE_CACHE && ret->size <= NBW_MAX_CACHED_SIZE);
    
    ret->alphabet = one->alphabet;
    ret->alphabet_size = one->alphabet_size;
    
    ret->char_labels = std::vector<std::string>(one->char_labels);
    
    // preserve pretty state labels if both input automata had them
    if(one->state_labels.size() >= one->size 
                            && two->state_labels.size() >= two->size){
        ret->state_labels = std::vector<std::string>();
        for(int i = 0; i < ret->size; i++)
            ret->state_labels.push_back(one->state_labels[i / two->size]
                                        + " & " 
                                        + two->state_labels[i % two->size]);
    }
    
    /* We allocate memory in which to store a bitvector of projected tracks --
     * although no tracks in this automaton should be projected yet!
     */
    ret->projected_tracks.resize(ret->alphabet_size);    

    /* A state in the product automaton is initial(final) if the corresponding
     * states in one and two are initial(final).
     */
    ret->initial.resize(ret->size);
    ret->final.resize(ret->size);
    for(int i = 0; i < ret->size; i++){
        if(one->initial[i / two->size] && two->initial[i % two->size])
            ret->initial.set(i);
        if(one->final[i / two->size] && two->final[i % two->size])
            ret->final.set(i);            
    }

    /* Set number of transitions */
    ret->num_transitions = one->num_transitions * two->size 
                            + two->num_transitions * one->size;
    
    /* allocate memory for transition matrix */
    ret->transition_matrix = new state_set_t[ (ret->alphabet_size) * (ret->size)];
    for(int i = 0; i < (ret->alphabet_size)*(ret->size); i++){
        ret->transition_matrix[i].resize(ret->size);
    }

    /* import transitions from input automata */
    for(int state = 0; state < ret->size; state++){
        for(int c = 0; c < ret->alphabet_size; c++){
            int s1 = state / two->size, s2 = state % two->size;
            state_set_t t1 = one->transition_matrix[s1*(one->alphabet_size) + c];
            state_set_t t2 = two->transition_matrix[s2*(two->alphabet_size) + c];
            for(int target = 0; target < ret->size; target++){
                int s1target = target / two->size;
                int s2target = target % two->size;
                if(t1[s1target] && t2[s2target])
                    ret->transition_matrix[state* (ret->alphabet_size) + c].set(target);        
            }
        }
    }

    /* Build a cache, if necessary.
     * (I don't expect this will be used frequently after product construction)
     */
    if(ret->use_cache){
        int cache_size = (ret->alphabet_size) * (1<<(ret->size));
        // allocate memory for transition cache
        ret->transition_cache = new state_set_t[cache_size];
        /* Build the entire cache.
         */
        for(unsigned long states = 0; states < (1 << (ret->size)); states++){
            state_set_t states_from(ret->size, states);
            for(int c = 0; c < ret->alphabet_size; c++){   
                int cache_index = (states*(ret->alphabet_size)) + c;
                ret->transition_cache[cache_index].resize( ret->size );
                for(int i = 0; i < ret->size; i++){
                    if(states_from[i])
                        ret->transition_cache[cache_index] |= ret->transition_matrix[ i * ret->alphabet_size + c];
                }
            }
        }
        
    }
    return ret;    
} // end NBW* NBW::product(NBW* one, NBW* two)


/** Returns true IFF the language of the automaton is empty. 
 */
bool NBW::is_empty(){
    this->trim();
    if(this->size > 1)
        return false;
    else{
        return !(this->initial[0] && this->final[0] && this->num_transitions > 0);
    }
} // end bool NBW::is_empty() const


void NBW::project(int track_index){
    this->trimmed = false;
    for(int c1 = 0; c1 < this->alphabet_size; c1++){
        int c2 = c1 ^ (1 << track_index); // toggle the bit corresponding to that track
        if(c2 > c1){
            for(int s = 0; s < this->size; s++){
                transition_matrix[s*alphabet_size + c1] |= transition_matrix[s*alphabet_size + c2];
                transition_matrix[s*alphabet_size + c2] |= transition_matrix[s*alphabet_size + c1];
            }
        }
    }
    
    if(this->use_cache){
        /* Rebuild the cache, if necessary.
         */
        for(unsigned long states = 0; states < (1 << (this->size)); states++){
            state_set_t states_from(this->size, states);
            for(int c = 0; c < this->alphabet_size; c++){   
                int cache_index = (states*(this->alphabet_size)) + c;
                for(int i = 0; i < this->size; i++){
                    if(states_from[i])
                        this->transition_cache[cache_index] |= this->transition_matrix[ i * this->alphabet_size + c];
                }
            }
        }
        
    }
} // end void NBW::project(int track_index)

state_set_t NBW::accessible_states() const{
    state_set_t accessible(this->initial);
    std::vector<int> bfs_queue;
    
    for(int i = 0; i < this->size; i++){
        if(this->initial[i])
            bfs_queue.push_back(i);
    }  
    
    for(int i = 0; i < bfs_queue.size(); i++){
        for(int c = 0; c < this->alphabet_size; c++){
            state_set_t targets = this->transition_matrix[bfs_queue[i] * this->alphabet_size + c];
            for(int new_state = 0; new_state < this->size; new_state++){
                if(targets[new_state] && !accessible[new_state]){
                    bfs_queue.push_back(new_state);
                    accessible.set(new_state);                
                }          
            }
        }
    }
    
    return accessible;
}

/* A coaccessible state is one with a path to a loop containing an accept state.
 */
state_set_t NBW::coaccessible_states() const{
    
    /* First step: identify accept states which are in loops. We'll call these
     * states "alive".
     */   
    state_set_t alive(this->size);

    // look for final states with loops (trivial SCCs)
    for(int i = 0; i < this->size; i++){
        if(this->final[i]){
            for(int j = 0; j < alphabet_size && !alive[i]; j++) {
                if(this->transition_matrix[i*alphabet_size + j][i]){
                    alive.set(i); // accept state in a self-loop is alive
                }
            }
        }        
    }         
    
    // build SCC's for the automaton, ignoring transitions
    BoostGraph g(this->size);
    for(int i = 0; i < this->size; i++){
        for(int j = 0; j < this->alphabet_size; j++){
            assert(this->transition_matrix[i*alphabet_size + j].size() == this->size);
            for(int target = 0; target < this->size; target++){
                if(this->transition_matrix[i*alphabet_size + j][target])
                    boost::add_edge(i, target, g);
            }
        }
    }
    
    std::vector<int> sccs(boost::num_vertices(g));
    int num_sccs = boost::strong_components(g, &sccs[0]);
    bool nontrivial_component[num_sccs]; //component already had size >= 1
    bool component_satisfied[num_sccs];  //already saw a final state in component
    
    for(int c = 0; c < num_sccs; c++){
        nontrivial_component[c] = false;
        component_satisfied[c] = false;
    }
    
    // Determine which of the connected components are alive.
    for(int i = 0; i < this->size; i++){
        if(this->final[i]){
            if(nontrivial_component[sccs[i]]){
                alive[i] = true;
            }
            component_satisfied[sccs[i]] = true;
        } else {
            // we already saw a final state in this component... we just needed
            // to see if there were any other states in the component.
            // there are, so the automaton is empty
            if(component_satisfied[sccs[i]]){
                alive[i] = true;
            }
        }
        nontrivial_component[sccs[i]] = true;        
    }

    /* At this point, the states that need to be marked as live are those
     * from which a live state can be reached, but are not in connected 
     * components (or were the first states seen in a connected component
     * containining a final state). 
     * We will run BFS in the reverse graph from the alive states
     * to find them.
     * Current algorithm (no BGL): reversing the graph: O(V^2).
     * BFS: O(V + E).
     */

    // reverse graph -- reverse_accessible[i] contains all the states which
    // transition to i.
    std::vector<int> reverse_accessible[this->size];
    for(int i = 0; i < this->size; i++){
        reverse_accessible[i] = std::vector<int>();
    }
    
    for(int j = 0; j < this->size; j++){
        for(int i = 0; i < this->size; i++){
            bool hit = false;
            for(int k = 0; k < this->alphabet_size && !hit; k++){            
                if(this->transition_matrix[j*this->alphabet_size + k][i]){
                    reverse_accessible[i].push_back(j);
                    hit = true;
                }
            }
        }
    }
    
    std::vector<int> search_queue;
    for(int i = 0; i < this->size; i++){
        if(alive[i])
            search_queue.push_back(i);
    }
    
    for(int i = 0; i < search_queue.size(); i++){
        for(int j = 0; j < reverse_accessible[search_queue[i]].size(); j++){
            if(!alive[reverse_accessible[search_queue[i]][j]]){
                alive.set(reverse_accessible[search_queue[i]][j]);
                search_queue.push_back(reverse_accessible[search_queue[i]][j]);
            }
        }
    }
    
    return alive;
}

/* Remove any states which are not (accessible and coaccessible) and
 * condense the automaton accordingly. Has no effect if the size of the
 * automaton is 1, so that empty machines can be represented by a single state.
 * Returns the number of states saved by this process (old size - new size).
 */
int NBW::trim(){
    if(this->trimmed)
        return 0;    
    
    if(this->size <= 1)
        return 0;

    // calculate states to keep
    state_set_t acc = this->accessible_states();
    state_set_t coacc = this->coaccessible_states();
    state_set_t keep = acc & coacc;
    
    // calculate the old labels for the states we're keeping
    int new_size = keep.count();
    
    /* If the machine is empty, we'll quicken this process by creating a
     * single empty state.
     */
    if(new_size == 0){
        int states_saved = this->size - 1;
        // Free all of the useless data structures.
        this->state_labels.clear(); // don't need this anymore!
        delete [] this->transition_matrix;
        if (this->use_cache)
            delete [] this->transition_cache;
        SafraTree::reset();
        this->size = 1;
        this->transition_matrix = new state_set_t[this->alphabet_size];
        for(int i = 0; i < this->alphabet_size; i++){
            this->transition_matrix[i].resize(1);
        }
        this->initial.resize(1);
        this->initial.set(0);
        this->final.resize(1);
        this->final.reset(0);
        this->num_transitions = 0;
    
        // I really hate caching. Deal with it if necessary...
        // Just for consistency's sake at least. Waste of memory though.
        this->use_cache = (NBW_USE_CACHE && this->size <= NBW_MAX_CACHED_SIZE);
        if(this->use_cache){
            int cache_size = (this->alphabet_size) * (2);
            // allocate memory for new transition cache
            this->transition_cache = new state_set_t[cache_size];
            for(int i = 0; i < cache_size; i++){
                this->transition_cache[i].resize(1);
            }
        } 
    
        return states_saved;
    }
    
    
    
    int old_labels[new_size];
    int new_label = 0;
    for(int i = 0; i < this->size; i++){
        if(keep[i]){
            old_labels[new_label] = i;
            new_label++;
        }
    }
    
    // allocate memory for a new transition matrix
    state_set_t* new_tm = new state_set_t[ (this->alphabet_size) * (new_size)];
    for(int i = 0; i < (this->alphabet_size)*(new_size); i++){
        new_tm[i].resize(new_size);
    }
    
    // Calculate the new transition matrix
    int new_transitions = 0;
    for(int s1 = 0; s1 < new_size; s1++){
        int old_state = old_labels[s1];
        for(int c = 0; c < this->alphabet_size; c++){
            for(int s2 = 0; s2 < new_size; s2++){
                if(this->transition_matrix[old_state*this->alphabet_size + c][old_labels[s2]]){
                    new_tm[s1*this->alphabet_size + c].set(s2);
                    new_transitions++;
                }
            }
        }
    }

    state_set_t new_initial(new_size);
    state_set_t new_final(new_size);
    for(int i = 0; i < new_size; i++){
        if(this->initial[old_labels[i]])
            new_initial.set(i);
        if(this->final[old_labels[i]])
            new_final.set(i);    
    }

    /* Free the old memory reserved for this automaton; we'll do it here since
     * this is the very first point at which it's no longer needed. This also
     * represents the point at which we start borking fields of the object --
     * transitions should be considered invalid from here to right before we 
     * return, and trimmed is set to true.
     * IMPORTANT / README / NOTE: Something needs to be (invented and then) 
     * locked at this point if you ever want to use this object in a 
     * multithreaded / multiple-reentrant context.
     */
    delete [] this->transition_matrix;
    if (this->use_cache)
        delete [] this->transition_cache;
    SafraTree::reset();

    /* Update the list of state labels if necessary. */
    if(this->state_labels.size() > 0){
        std::vector<std::string> new_labels;
        for(int i = 0; i < this->size; i++){
            if(keep[i])
                new_labels.push_back(state_labels[i]);
        }
        this->state_labels.clear();
        this->state_labels = new_labels;
    }
    
    int states_saved = this->size - new_size;

    this->size = new_size;
    this->transition_matrix = new_tm;
    this->num_transitions = new_transitions;
    this->initial = new_initial;
    this->final = new_final;

    // Deal with caching if necessary.
    this->use_cache = (NBW_USE_CACHE && this->size <= NBW_MAX_CACHED_SIZE);
    if(this->use_cache){
        int cache_size = (this->alphabet_size) * (1<<(this->size));
        // allocate memory for new transition cache
        this->transition_cache = new state_set_t[cache_size];
        for(unsigned long states = 0; states < (1 << (this->size)); states++){
            state_set_t states_from(this->size, states);
            for(int c = 0; c < this->alphabet_size; c++){   
                int cache_index = (states*(this->alphabet_size)) + c;
                this->transition_cache[cache_index].resize( this->size );
                for(int i = 0; i < this->size; i++){
                    if(states_from[i])
                        this->transition_cache[cache_index] |= this->transition_matrix[ i * this->alphabet_size + c];
                }
            }
        } 
    }    

    this->trimmed = true;

    return states_saved;
}


/** DRW.cpp: implements a deterministic Rabin automaton which
 *  recognizes languages of (one-way) infinite words.
 *
 *  (c) Joe Gershenson, 2009
 */
 
// std libraries for strings, I/O, and vector.
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>


#include <boost/graph/strong_components.hpp>
#include <boost/graph/transitive_closure.hpp>

#include "DRW.hpp"
#include "utils.hpp"

std::string DRW::to_string(){
    using namespace std;
    string s;
    s.reserve(300);    
    s.append("#----- omega-automaton (DRW) ----- \n");
    s.append("RABIN\n");
    s.append("# Number of states: \n");
    s.append(INT_TO_STR(this->size));
    s.append("\n# Size of alphabet: \n");
    s.append(INT_TO_STR(this->alphabet_size));
    s.append("\n# List of transitions: \n");
    for(int i = 0; i < this->size; i++){
        for(int j = 0; j < this->alphabet_size; j++){
            int k = this->transition_matrix[i][j];
            s.append( INT_TO_STR(i+1) );
            s.append( " > " );
            s.append( INT_TO_STR(j+1) );
            s.append( " > " );
            s.append( INT_TO_STR(k+1) );
            s.append( "\n" );
    	}
    }
    s.append("# Initial state\n");
    s.append(INT_TO_STR(this->initial_state + 1));
    s.append("\n");
    s.append("# Rabin pairs: \n");
    for(int i = 0; i < this->pairs.size(); i++){
        for(int j = 0; j < this->size; j++){
            if(this->pairs[i]->finite[j]){
                s.append( INT_TO_STR(j+1) );
                s.append(" ");
            }
        }
        s.append("| ");
        for(int j = 0; j < this->size; j++){
            if(this->pairs[i]->infinite[j]){
                s.append( INT_TO_STR(j+1) );
                s.append(" ");
            }
        }
        s.append("\n");    
    }
    s.append("# EOF\n");
    return s;
}


std::string DRW::to_GASt_string(){
    using namespace std;
    if(!SAVE_TREE_DATA)
        return string("Tree data not saved -- set SAVE_TREE_DATA to true in SafraTree.hpp to use this feature");
    string s;
    s.append("Deterministic Rabin-Automaton according to Safra:\n");
    s.append("\n");
    s.append(INT_TO_STR(this->size));
    s.append(" states:\n");
    for(int i = 0; i < this->size; i++){
        s.append("s");
        s.append(INT_TO_STR(i+1));
        s.append(":\n");
        s.append(SafraTree::get_tree(i)->to_GASt_string());
    }
    return s;
}

	/*
	 * Return a visual representation of the transition graph suitable
	 * for processing with dot (graphviz.org).
	 */
std::string DRW::to_digraph(){
    using namespace std;

    string s;
    s.append("digraph rabin_automaton {\n");
    s.append("    node [shape=circle];\n");

    // initial state gets an arrow in from an invisible state
    s.append("    initial_invis [style=invis];\n");
    s.append("    initial_invis ->");
    s.append(INT_TO_STR(this->initial_state + 1));
    s.append(";\n");
	
    for(int i = 0; i < this->size; i++){
        for(int k = 0; k < this->size; k++){
            string label;
            // write down all characters on which state i goes to state k
            for (int j = 0; j < this->alphabet_size; j++){
                int target = this->transition_matrix[i][j];
                if (target == k){
                    label.append(",");
                    if(this->char_labels.size() == 0)
                        label.append(1, this->alphabet[j]);
                    else
                        label.append(this->char_labels[j]);
                }
            }
            if(label.size() > 0){
                s.append("    ");
                s.append( INT_TO_STR(i+1) );
                s.append(" -> ");
                s.append( INT_TO_STR(k+1) );
                s.append(" [label=\"");
                s.append(label.erase(0,1));
                s.append("\"];\n");
            }
        }
    }
    s.append("}");
    return s;
}

int DRW::transition(int state, int character) const{
    return this->transition_matrix[state-1][character-1];
}

DRW::DRW(){
    this->alphabet = default_alphabet; //defined in utils.cpp
    this->char_labels.resize(0);
    this->t_matrix = NULL;
    this->sccs = NULL;
    /* The field num_sccs is used to track whether the boost components
     * have been initialized or not.
     */
    this->num_sccs = -1; 
}

DRW::~DRW(){
    delete this->t_matrix;
    delete this->sccs;
    for(int i = 0; i < this->transition_matrix.size(); i++)
        delete [] this->transition_matrix[i];
    for(int i = 0; i < this->pairs.size(); i++)
        delete this->pairs[i];
}

RabinPair::RabinPair(int size){
    this->infinite.resize(size);
    this->finite.resize(size);
}

RabinPair::~RabinPair(){
    this->infinite.clear();
    this->finite.clear();
}

DRW* DRW::parse(char* filename) {
    using namespace std;
    string s;

    ifstream inf (filename);
    if (inf.good()) {
        get_next_line(inf, s);


        // check to see if file is in GASt format
        if( (s.compare("RABIN") != 0) ){
            cout << "IO Error! Check file format." << endl;
            return NULL;
        }
        
        // allocate memory
        DRW* ret = new DRW();

        // determine automaton size (number of states)
        get_next_line(inf, s);
        ret->size = atoi(s.c_str());
        
        // determine alphabet size
        get_next_line(inf, s);
        ret->alphabet_size = atoi(s.c_str());
        
        // allocate memory for transition matrix
        for(int i = 0; i < ret->size; i++){
            ret->transition_matrix.push_back(new int[ret->alphabet_size]);
        }        

        // process transitions
        for(int i = 0; i < (ret->size * ret->alphabet_size); i++){
            get_next_line(inf, s);
            string part;
            istringstream iss(s);
            getline(iss, part, ' ');
            int state_from = atoi(part.c_str());
            getline(iss, part, ' '); // now part contains the arrow >
            getline(iss, part, ' '); // part now contains the character
            int char_on = atoi(part.c_str());
            getline(iss, part, ' '); // part contains the second arrow
            getline(iss, part, ' '); // part contains the target state
            int state_to = atoi(part.c_str());
            ret->transition_matrix[state_from-1][char_on-1] = state_to-1;
        }
        
        // read initial state
        get_next_line(inf, s);        
        ret->initial_state = atoi(s.c_str())-1;
        
        
        // allocate memory for final state set
        get_next_line(inf, s);
        
        while(s.size() > 0){
            // Read a Rabin pair
            RabinPair* rp = new RabinPair(ret->size);
            string part;
            istringstream iss(s);
            getline(iss, part, ' ');
            // read the finitely-many-times states
            while(part.compare("|")){
                rp->finite.set(atoi(part.c_str())-1);
                getline(iss, part, ' ');
            }
            // throw away the '|' in the middle
            getline(iss, part, ' ');
            // read the infinitely-many-times states            
            while(part.size() > 0){
                rp->infinite.set(atoi(part.c_str())-1);
                getline(iss, part, ' ');
            }
            // add this RAbin pair to the automaton
            ret->pairs.push_back(rp);
        
            // get ready to read the next Rabin pair
            get_next_line(inf, s);
        }
        std::vector<RabinPair*> pairs;

        
        return ret;
    }
    // code reaching this point implies IO error.
    std::cout << "I/O error with " << filename << "; aborting.\n";
    return 0;
}

/*
 * Build a version of the transition matrix using the 
 * Boost Graph Library (http://www.boost.org), and calculate the
 * strongly connected components of the graph. Information is stored in the 
 * private fields @field t_matrix (for the graph itself) and @field sccs
 * (for the strongly connected components).
 */
void DRW::build_boost_components(){
    using namespace boost;
    
    if(this->num_sccs == -1){ 
        // delete any old data as necessary
        if(this->t_matrix != NULL)
            delete this->t_matrix;
        if(this->sccs != NULL)
            delete this->sccs;    
            
        this->t_matrix = new BoostGraph(this->size);
        
        for(int i = 0; i < this->size; i++){
            for(int j = 0; j < this->alphabet_size; j++){
                add_edge(i, this->transition_matrix[i][j], *t_matrix);
            }
        }
        
        this->sccs = new std::vector<int>(num_vertices(*(this->t_matrix)));
        this->num_sccs = strong_components(*(this->t_matrix), &(*(this->sccs))[0]);    
    }
}

/* 
 * Test function: print all strongly connected components of the graph.
 */
void DRW::print_components(){
  this->build_boost_components();
  using namespace boost;   
  std::cout << "Total number of components: " << this->num_sccs << std::endl;
  std::vector<int>::size_type i;
  for (i = 0; i != this->sccs->size(); ++i)
    std::cout << "State " << (i+1) <<" is in component " << (*(this->sccs))[i] << std::endl;
}
        
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
bool DRW::is_empty(){
    using namespace boost;  
    
    for(int p = 0; p < this->pairs.size(); p++){
        /* To check the satifiability of the current pair, we first check to see
         * if there is a state in INF with a self-loop. This satisfies the pair
         * if present. Otherwise, we build a copy of 
         * the reachability graph (transition matrix ignoring labels) for the 
         * automaton, ignoring the states in FIN for the current pair. Then,
         * we calculate the strongly connected components of that graph. 
         * If there is any component of size >1 with a node in INF, the pair is
         * satisfied.
         */
        
        // look for single states with loops
        for(int i = 0; i < this->size; i++){
            if(this->pairs[p]->infinite[i]){
                for(int j = 0; j < this->alphabet_size; j++)
                    if(this->transition_matrix[i][j] = i)
                        return false; //automaton is not empty
            }        
        }         
        
        // build SCC's for (G \ FIN) 
        BoostGraph g(this->size);
        for(int i = 0; i < this->size; i++){
            if(!this->pairs[p]->finite[i]){
                for(int j = 0; j < this->alphabet_size; j++){
                    int target = this->transition_matrix[i][j];
                    if(!this->pairs[p]->finite[target])
                        add_edge(i, target, g);
                }
            }
        }
        
        std::vector<int> sccs(num_vertices(g));
        int num_sccs = strong_components(g, &sccs[0]);
        bool nontrivial_component[num_sccs];
        
        for(int c = 0; c < num_sccs; c++)
            nontrivial_component[c] = false;
        
        for(int i = 0; i < this->size; i++){
            if(this->pairs[p]->infinite[i] && nontrivial_component[sccs[i]])
                return false;
            nontrivial_component[sccs[i]] = true;        
        }
    }    
    
    return true;
}


/*
 *  TODO: still completely borked.
 * Determine if the language of the automaton is universal.

 */
bool DRW::is_universal(){
    return false;
}
        
/* 
 * Generate and return a BŸchi automaton which accepts the complement 
 * of the language accepted by this automaton.
 */
NBW* DRW::complement(){

    std::vector< boost::tuple<int, int, int> > adjacency_list;
    std::vector< CompState* > seen; // keep track of states we've seen
    std::vector< CompState* > work_queue;

    // add initial state
    // note that states in the initial part (p,0,0) need no memory for statesets
    CompState* initial = new CompState(0);
    initial->rabin_state = this->initial_state;
    initial->in_initial_part = true;
    initial->buchi_index = 0;
    work_queue.push_back(initial);
    seen.push_back(initial);
        
    // calculate reachable part of automaton
    while(!work_queue.empty()){
        std::vector<CompState*> new_work_queue;
        int max = work_queue.size();
        
        for(int i = 0; i < max; i++){
            CompState* current = work_queue[i];
            for(int a = 0; a < this->alphabet_size; a++){
                int q = this->transition_matrix[work_queue[i]->rabin_state][a];
                if(current->in_initial_part){
                    // (p,0,0) -a-> (q,0,0)
                    CompState* next_state = new CompState(0);
                    next_state->rabin_state = q;
                    next_state->in_initial_part = true;
                    int index = next_state->get_index(seen);
                    if(index == -1){
                        new_work_queue.push_back(next_state);
                        index = next_state->buchi_index;
                    } else {
                        delete next_state; // clean up extra object
                        next_state = seen[index];
                    }
                    adjacency_list.push_back(boost::make_tuple(work_queue[i]->buchi_index, a, index));
                    
                    // (p,0,0) -a-> (q,\0,\0)
                    next_state = new CompState(this->pairs.size());
                    next_state->rabin_state = q;
                    next_state->in_initial_part = false;
                    index = next_state->get_index(seen);
                    if(index == -1){
                        new_work_queue.push_back(next_state);
                        index = next_state->buchi_index;
                    } else {
                        delete next_state; // clean up extra object
                        next_state = seen[index];
                    }
                    adjacency_list.push_back(boost::make_tuple(work_queue[i]->buchi_index, a, index));                    
                
                } else {
                    // (p,s1,s2) -a-> (q,s1',s2')
                    // where s1 tracks finite pairs hit, s2 infinite pairs hit
                    // and s1 hits cancel s2 hits (leaving the pair unsatisfied).
                    // Since a state is final if s2 is empty, the complement
                    // machine must keep all pairs unsatisfied forever
                    CompState* next_state = new CompState(this->pairs.size());
                    next_state->rabin_state = q;
                    next_state->in_initial_part = false;
                    next_state->s1 = state_set_t( current->s1 );
                    next_state->s2 = state_set_t( current->s2 );
                    
                    for(int pair = 0; pair < this->pairs.size(); pair++){
                        if(this->pairs[pair]->finite[q])
                            next_state->s1.set(pair);
                        else if(this->pairs[pair]->infinite[q])
                            next_state->s2.set(pair);
                    }
                    
                    if(next_state->s2.is_subset_of(next_state->s1)){
                        next_state->s1 -= next_state->s2;
                        next_state->s2.reset();   
                    }
                    
                    int index = next_state->get_index(seen);                    
                    
                    if(index == -1){
                        new_work_queue.push_back(next_state);
                        index = next_state->buchi_index;
                    } else {
                        delete next_state; // clean up extra object
                        next_state = seen[index];
                    }               

                    adjacency_list.push_back(boost::make_tuple(work_queue[i]->buchi_index, a, next_state->buchi_index));                                        

                
                }    
            } // end for each character in the alphabet
        } // end for each state in the work queue
        work_queue = new_work_queue;
        
    } // end calculating reachable part of new automaton
    
    
    /* We have finished calculating the adjacency list; 
       now calculate the rest of the data to build the automaton. */    
        
    int nbw_size = seen.size();
    int nbw_alphabet_size = this->alphabet_size;
    int num_transitions = adjacency_list.size();
    
    // copy character semantics table
    std::vector<std::string> nbw_char_labels(this->char_labels);
    
    // Create a state semantics table (state labels) corresponding
    // to the complementation procedure we used.
    std::vector<std::string> nbw_state_labels;
    
    for(int i = 0; i < nbw_size; i++){        
        if(seen[i]->in_initial_part){
            std::string label("(");
            label += INT_TO_STR(seen[i]->rabin_state + 1);
            label += ", initial)";
            nbw_state_labels.push_back(label);
        } else {
            std::string label("(");
            label += INT_TO_STR(seen[i]->rabin_state + 1);
            label += ",";
            std::string s1; std::string s2;
            boost::to_string(seen[i]->s1, s1);
            boost::to_string(seen[i]->s2, s2);
            label += (s1 + "," + s2 + ")");
            nbw_state_labels.push_back(label);    
        }
    } 
    
    // set initial and final states.
    // Note that the initial state is state 0 since we added it first.
    state_set_t nbw_initial(nbw_size);
    nbw_initial.set(0);

    state_set_t nbw_final(nbw_size);
    for(int i = 0; i < nbw_size; i++){ //(ret->size = seen.size())
        if(!seen[i]->in_initial_part){
            if(seen[i]->s2.none()){
                nbw_final.set(i);
            }
        }
    }
    
    // Delete CompStates to free memory
    for(int i = 0; i < seen.size(); i++)
        delete seen[i];
    
    NBW* ret = new NBW(nbw_size, 
                   nbw_alphabet_size, 
                   adjacency_list,
                   nbw_initial,
                   nbw_final,
                   nbw_char_labels,
                   nbw_state_labels);
    
    return ret;
}


/***************** Implementation of private class CompState ************/
/* Used in DRW::complement to build states of a Buchi automaton accepting
 * the complement of the language
 */

DRW::CompState::CompState(int state_set_size){
    this->s1.resize(state_set_size);
    this->s2.resize(state_set_size);
}

/** Simple equality check -- if all the fields are equal they are equal.
 */
bool DRW::CompState::operator==(const CompState& other){
    return (this->rabin_state == other.rabin_state) 
        && (this->in_initial_part == other.in_initial_part)
        && (this->in_initial_part || ((this->s1 == other.s1) && (this->s2 == other.s2)));
}

/** Get the index of an equivalent CompState in the vector, 
 *  if it is in the vector at all.
 *  Returning -1 indicates that we just added it... you should retrieve the new
 *  reference. Returning any other value indicates that it is a redundant object...
 *  you should delete this object.
 */
int DRW::CompState::get_index(std::vector<CompState*>& seen){
    for(int i = 0; i < seen.size(); i++)
        if(*(seen[i]) == *(this))        
            return i;
    this->buchi_index = seen.size();
    seen.push_back(this);
    return -1;
}
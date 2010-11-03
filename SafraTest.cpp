/** @file SafraTest.cpp
 *  Code for testing Safra's construction.
 *
 *  (c) Joe Gershenson, 2009
 */

#include <string>
#include <iostream>
#include <fstream>

#include <boost/unordered_map.hpp>

#include <omp.h>
#include <time.h>

#include <stdio.h>
#include <stdlib.h>

#include "SafraTest.hpp"
#include "SafraTree.hpp"
#include "NBW.hpp"

using namespace std;

int main(int argc, char** argv){
    determinize(argv[1], argv[2]);
    exit(0);
}

void determinize(char* infile, char* outfile){
    NBW* nbw = NBW::parse(infile);
    time_t start = time(NULL);
    DRW* drw = nbw->determinize();
    time_t end = time(NULL);
    std::ofstream out;
    out.open(outfile);
    out << drw->to_string();
    out.flush();
    out.close();
    std::cout << "  Rabin states: " << (drw->size) << std::endl;
    std::cout << " PLZPLZPLZ \n" << (drw->complement()->to_digraph()) << std::endl;
    delete nbw;
    delete drw;
}

int run_random_trials(int argc, char** argv){
    if(argc < 6){
        cout << "Usage: " << argv[0] << " automaton-size alphabet-size transition-density final-density num-trials [report-freq]" << endl;
    } else {
        int SIZE = atoi(argv[1]);
        int ALPHABET_SIZE = atoi(argv[2]);
        double TRANS_PROB = atof(argv[3]);
        double FINAL_PROB = atof(argv[4]);
    
        int NUM_TRIALS = atoi(argv[5]);
        
        bool report_progress;
        int REPORT_FREQ;
        if(argc < 7)
            report_progress = false;
        else {
            report_progress = true;
            REPORT_FREQ = atoi(argv[6]);
        }
    
        cout << "Determinizing random automata." << endl;
        // Seed RNG
        time_t seed = time(NULL);
        cout << " -- RNG seed is " << seed << endl;
        srand(seed);
        
        cout << " -- Size: " << SIZE << endl;
        cout << " -- Alphabet size: " << ALPHABET_SIZE << endl;
        cout << " -- Transition probability: " << TRANS_PROB << endl;
        cout << " -- Final state probability: " << FINAL_PROB << endl;

        cout << endl << endl;
        
        time_t dettime = 0;
        time_t gentime = 0;
        time_t freetime = 0;
        
        unsigned int max = 0;
        unsigned int min = 1<<31;
        unsigned int sum = 0;
    
        for(int trial = 0; trial < NUM_TRIALS; trial++){
            if(report_progress && (trial % REPORT_FREQ == 0))
                cout << " (" << ((100.0 * (trial)) / NUM_TRIALS) << "%) " << trial << " / " << NUM_TRIALS << endl;
        
            gentime -= time(NULL);
            NBW::NBW* my_nbw = NBW::build_random_automaton(SIZE, ALPHABET_SIZE, TRANS_PROB, FINAL_PROB);
            gentime += time(NULL);
            dettime -= time(NULL);
            DRW::DRW* drw = my_nbw->determinize();                    
            dettime += time(NULL);
            int size = drw->size;
            
            if(size > max)
                max = size;
            if(size < min)
                min = size;
            sum += size;
            
            freetime -= time(NULL);
            delete my_nbw;
            delete drw;
            freetime += time(NULL);
        } 
        
        cout << endl << endl;
        cout << "Results (" << NUM_TRIALS << " trials) :" << endl;
        cout << "  Max size: " << max << endl;
        cout << "  Min size: " << min << endl;
        cout << "  Avg size: " << (sum / NUM_TRIALS) << endl;
        cout << "  Random generation time: " << gentime << endl;
        cout << "  Determinization time: " << dettime << endl;
        cout << "  Deallocation time: " << freetime << endl;
        
        return 0;
    }
}

int print_averages(int argc, char** argv){

    int SIZE = atoi(argv[1]);
    int ALPHABET_SIZE = atoi(argv[2]);
    double TRANS_PROB = atof(argv[3]);
    double FINAL_PROB = atof(argv[4]);

    int NUM_TRIALS = atoi(argv[5]);
    

    // Seed RNG
    time_t seed = time(NULL);
    srand(seed);
    
    time_t dettime = 0;
    unsigned int sum = 0;

    for(int trial = 0; trial < NUM_TRIALS; trial++){
        NBW::NBW* my_nbw = NBW::build_random_automaton(SIZE, ALPHABET_SIZE, TRANS_PROB, FINAL_PROB);
        dettime -= time(NULL);
        DRW::DRW* drw = my_nbw->determinize();                    
        dettime += time(NULL);            
        sum += drw->size;            
        delete my_nbw;
        delete drw;
    } 
    
    cout << TRANS_PROB << " " << FINAL_PROB << " " << (sum / NUM_TRIALS) << " " << (dettime / NUM_TRIALS) << endl;
    
    return 0;
}

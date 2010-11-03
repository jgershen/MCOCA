/** @file gen_test.cpp
 *  Code for testing generation of Buechi automata to check CA properties.
 *
 *  (c) Joe Gershenson, 2009
 */

#include <string>
#include <iostream>
#include <fstream>

#include <stdio.h>
#include <stdlib.h>

#include <time.h>

#include "logic.hpp"
#include "buchi_gen.hpp"
#include "NBW.hpp"
#include "fol_parser.hpp"

using namespace std;

#define MIN_DEPTH 0
#define MAX_DEPTH 3

void graph_to_file(NBW* nbw, std::string outfile){
    std::ofstream out;
    out.open(outfile.c_str());
    out << nbw->to_digraph();
    out.flush();
    out.close();
}

void graph_to_file(DRW* nbw, std::string outfile){
    std::ofstream out;
    out.open(outfile.c_str());
    out << nbw->to_digraph();
    out.flush();
    out.close();
}


int check_nilpotency(int eca, bool debug, int max_level=MAX_DEPTH){
    //cout << "trying ECA " << INT_TO_STR(eca) << " level " << INT_TO_STR(level) << endl;
    for(int level = MIN_DEPTH; level < max_level; level++){
        Conjunction f;
        //cout << "trying ECA " << INT_TO_STR(eca) << " level " << INT_TO_STR(level) << endl;
        //cout << "whoo" << endl;
        
        for(int i = 0; i < level; i++){
            Literal* p = new Literal('a'+i, 'a'+i+1, eca, false);
            f.literals.push_back(p);
        }
        
        Literal* p = new Literal('a'+level, 'a'+level, eca, false); // last track is fixpoint 
        f.literals.push_back(p);       

        Quantifier forall(Q_POSITIVE,Q_UNIVERSAL,'a');
        f.add_inner_quantifier(forall);

        for(int i = 1; i <= level; i++){
            Quantifier q(Q_POSITIVE,Q_EXISTENTIAL,'a'+i);
            f.add_inner_quantifier(q);
        }

        if(debug)
            cout << "  " << f.to_string() << endl;
        
        std::vector<Conjunction*> form;
        form.push_back(&f);
        
        NBW* nbw = NBW::build_automaton(form, OMEGA);
        
        if(debug)
            graph_to_file(nbw, std::string("nbw.dot"));

        
        bool success = !(nbw->is_empty());

        delete nbw; //delete comp;

        if(success){
            cout << INT_TO_STR(eca) << " " << INT_TO_STR(level) << " (proved " << f.to_string() << " )" << endl;
            return level;
        }
    }
    return -1;
} // end quick_check_nilpotency(int eca)

bool check_surjectivity(int eca, bool debug){    
    // Ax Ey y->x
    // = ~ Ex ~ Ey y->x
    // = isEmpty (Ex ~ Ey y->x)


    Conjunction f;
        
    Literal* p1 = new Literal('a'+0, 'a'+1, eca, false);
    f.literals.push_back(p1);

    Quantifier u(Q_POSITIVE,Q_UNIVERSAL,'a'+1);
    f.add_inner_quantifier(u);

    Quantifier q(Q_POSITIVE,Q_EXISTENTIAL,'a'+0);
    f.add_inner_quantifier(q);

    std::vector<Conjunction*> form;
    form.push_back(&f); 


    NBW* nbw = NBW::build_automaton(form, OMEGA);


    if(!(nbw->is_empty())){
        cout << INT_TO_STR(eca) << " (proved " << f.to_string() << " )" << endl;
        delete nbw;
        return true;
    } else {
        delete nbw;        
        return false;
    }
} // end check_surjectivity

bool check_injectivity(int eca, bool debug){    
    Conjunction f;
        
    Literal* p1 = new Literal('a'+2, 'a'+0, eca, false);
    f.literals.push_back(p1);
    
    Literal* p2 = new Literal('a'+1, 'a'+0, eca, false);
    f.literals.push_back(p2);
    
    Literal* np2 = new Literal('a'+1,'a'+2, IDENTITY_ECA_NUM, true);
    f.neg_literals.push_back(np2); 
        
    Quantifier outer(false,false,'a'+0);
    f.add_inner_quantifier(outer);

    Quantifier middle(false,false,'a'+1);
    f.add_inner_quantifier(middle);
    
    Quantifier inner(false,false,'a'+2);
    f.add_inner_quantifier(inner);
    
    std::vector<Conjunction*> form;
    form.push_back(&f); 
 
    NBW* nbw = NBW::build_automaton(form, OMEGA);
    
    if(debug){
        string filename = "/Users/jgershen/research/r2/junk/nbw.dot";
        graph_to_file(nbw, filename);        
        cout << nbw->to_string() << endl;    
    }
            
    if(nbw->is_empty()){
        cout << INT_TO_STR(eca) << " (disproved " << f.to_string() << " )" << endl;

        delete nbw;
        
        return true;
    } else {
        
        delete nbw;
        
        return false;
    }
    
}

int main(int argc, char** argv){

    int nilpotent[256];
    bool injective[256];
    bool surjective[256];

    cout << "Nilpotent ECAs (with nilpotency index):" << endl;
    for(int eca = 0; eca < 256; eca++){
        nilpotent[eca] = check_nilpotency(eca, false);
    }
    
    cout << "Injective ECAs:" << endl;
    for(int eca = 0; eca < 256; eca++){
        injective[eca] = check_injectivity(eca, false);
    }

    
    cout << "Surjective ECAs:" << endl;
    for(int eca = 0; eca < 256; eca++){
        surjective[eca] = check_surjectivity(eca, false);
    }
    
    cout << "Nilpotency\n";
    cout << "{";
    for(int i = 0; i < 256; i++){
        cout << INT_TO_STR(nilpotent[i]);
        if(i < 255) cout << ",";
    }
    cout << "}\n";

    cout << "Injectivity\n";
    cout << "{";
    for(int i = 0; i < 256; i++){
        if (injective[i] )
            cout << "1" ;
        else
            cout << "0";
        if(i < 255) cout << ",";
    }
    cout << "}\n";
    
    cout << "Surjectivity\n";
    cout << "{";
    for(int i = 0; i < 256; i++){
        if (surjective[i])
            cout << "1" ;
        else
            cout << "0";
        if(i < 255) cout << ",";
    }
    cout << "}\n";
    
}



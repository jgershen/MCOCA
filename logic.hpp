#pragma once
#ifndef _LOGIC_H_
#define _LOGIC_H_

#include <string>
#include <vector>

#include <boost/tuple/tuple.hpp>
#include <boost/dynamic_bitset.hpp>

#include "utils.hpp"

#define Q_NEGATED true
#define Q_POSITIVE false
#define Q_EXISTENTIAL false
#define Q_UNIVERSAL true

class SymbolTable
{
    private:
        static std::vector<char> symbol_table;
    public:
        static int var_count();
        static int index_of(char var);  
        static char lookup(int index);
};

/* A Literal is an assertion that track at index i1 
 * goes to index i2 under a given ECA
 */
class Literal
{
    private:
        static int default_eca;
    public:
        bool negated;

        int i1;  // track 1 index
        int i2;  // track 2 index
        boost::dynamic_bitset<unsigned long> eca; // use 204 for equality
   
        void set_eca(int eca_num);
   
        Literal();
        Literal(char v1, char v2);
        Literal(char v1, char v2, int eca_num, bool negated);

        bool check(boost::dynamic_bitset<unsigned long> x, 
                   boost::dynamic_bitset<unsigned long> y, 
                   boost::dynamic_bitset<unsigned long> z);  
                   
        std::string to_string() const;
        
        static void set_default_eca(int eca_num);
};

class Quantifier
{
    public:
        bool negated;
        bool universal;
        int variable_index;
        std::string to_string() const;        
        
        Quantifier(bool negated, bool universal, char var);
};

class Conjunction
{
    public:    

        std::vector<Literal*> literals;
        std::vector<Quantifier> quantifiers;
        std::vector<Literal*> neg_literals;
    
   
        static std::vector<Conjunction*>* last_formula_parsed;
    
        bool check(boost::dynamic_bitset<unsigned long> x, 
                   boost::dynamic_bitset<unsigned long> y, 
                   boost::dynamic_bitset<unsigned long> z);

        std::string to_string() const;
         
        void add_literal(Literal* l);
        void add_inner_quantifier(const Quantifier& q);
        void add_outer_quantifier(const Quantifier& q);
        
        Conjunction();
        //Conjunction(int max_tracks);
        //Conjunction(const std::vector<Quantifier>& global_quantifiers, int max_tracks);
};

#endif
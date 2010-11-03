#include "logic.hpp"
#include <iostream>

std::vector<Conjunction*>* Conjunction::last_formula_parsed;
std::vector<char> SymbolTable::symbol_table = std::vector<char>();

int Literal::default_eca = 0;

int SymbolTable::var_count(){
    return symbol_table.size();
}

int SymbolTable::index_of(char var){
    for(int i = 0; i < symbol_table.size(); i++){
        if(symbol_table[i] == var) {
            return i;                       
        }
    }
    symbol_table.push_back(var);
    return symbol_table.size() - 1;
}

char SymbolTable::lookup(int index){
    assert(index < symbol_table.size());
    return symbol_table[index];
}

Quantifier::Quantifier(bool negated, bool universal, char var){
    this->negated = negated;
    this->universal = universal;
    this->variable_index = SymbolTable::index_of(var);
}

std::string Quantifier::to_string() const{
    std::string value;
    if(this->negated)
        value += "~";
    if(this->universal)
        value += "A";
    else
        value += "E";
    value += SymbolTable::lookup(this->variable_index);
    return value;
}

/**************** Implementation of Conjunction *************************/

void Conjunction::add_literal(Literal* l){
    if(l->negated)
        this->neg_literals.push_back(l);
    else
        this->literals.push_back(l);
}

void Conjunction::add_inner_quantifier(const Quantifier& q){
    this->quantifiers.push_back(q);
}

void Conjunction::add_outer_quantifier(const Quantifier& q){
    this->quantifiers.insert(this->quantifiers.begin(), Quantifier(q));
}

/** Check to see if every positive literal in a conjunction is satisfied. Does
 *  not check negative literals.
 */
bool Conjunction::check(boost::dynamic_bitset<unsigned long> x, 
           boost::dynamic_bitset<unsigned long> y, 
           boost::dynamic_bitset<unsigned long> z){
    for(int i = 0; i < this->literals.size(); i++){
        if(! (this->literals[i]->check(x,y,z)) )
            return false;
    }
    return true;
}

std::string Conjunction::to_string() const{
    std::string value;
    for(int i = 0; i < this->quantifiers.size(); i++){
        value += quantifiers[i].to_string();
        value += " ";
    }

    value += ("(");
    for(int i=0; i < this->literals.size(); i++){
        value += this->literals[i]->to_string();
        if( (i+1) < this->literals.size() || this->neg_literals.size() > 0)
            value += " & ";        
    }
    for(int i=0; i < this->neg_literals.size(); i++){
        value += this->neg_literals[i]->to_string();
        if( (i+1) < this->neg_literals.size())
            value += " & ";        
    }
    
    value += std::string(")");
    return value;
}

Conjunction::Conjunction(){
    this->literals = std::vector<Literal*>();
    this->neg_literals = std::vector<Literal*>();    
}

/*************** Implementation of Literal ************************/

void Literal::set_default_eca(int eca_num){
    Literal::default_eca = eca_num;
}


/* Checks whether the condition expressed in this literal holds between the
 * specified slices. Ignores the negation status of this literal. Ex: if
 * this literal is "the tracks are equal" and the tracks *are* equal
 * this function will return true. If this literal is "NOT(the tracks are equal)"
 * and the tracks *are* equal, it will still return true.
 */
bool Literal::check(boost::dynamic_bitset<unsigned long> x, 
           boost::dynamic_bitset<unsigned long> y, 
           boost::dynamic_bitset<unsigned long> z){
    bool applies = eca[(4 * x[i1]) + (2 * y[i1]) + (z[i1])] == y[i2];
    return applies;
}

void Literal::set_eca(int eca_num){
    this->eca = boost::dynamic_bitset<unsigned long>(8, eca_num);
}

/*
 * Factory style constructor. If you call it, manually set fields before 
 * doing anything with the object!
 */
Literal::Literal(){
    this->i1 = -1;
    this->i2 = -1;
    this->negated = true;
}

Literal::Literal(char v1, char v2){
    int i1 = SymbolTable::index_of(v1);
    int i2 = SymbolTable::index_of(v2);
    this->i1 = i1;
    this->i2 = i2;
    this->eca = boost::dynamic_bitset<unsigned long>(8, default_eca);
    this->negated = false;
}

Literal::Literal(char v1, char v2, int eca_num, bool negated){
    int i1 = SymbolTable::index_of(v1);
    int i2 = SymbolTable::index_of(v2);
    this->i1 = i1;
    this->i2 = i2;
    this->eca = boost::dynamic_bitset<unsigned long>(8, eca_num);
    this->negated = negated;
}

std::string Literal::to_string() const{
    std::string ret;
    if(negated)
        ret += "~(";
    ret += SymbolTable::lookup(i1);
    if(this->eca.to_ulong() == IDENTITY_ECA_NUM)
        ret += "==";
    else
        ret += "->";
    
/*    ret += ">";
    ret += INT_TO_STR(this->eca.to_ulong());
    ret += ">";
*/
    ret += SymbolTable::lookup(i2);
    if(negated)
        ret += ")";
    return ret;
}

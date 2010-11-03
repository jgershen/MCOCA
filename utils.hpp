/**  @file utils.hpp
 *   @author Joe Gershenson
 *   
 *   Contains prototypes for utility functions that are used in multiple files 
 *   or classes. Also includes #define statements or extern declarations which 
 *   are used across multiple files. The goal is to concentrate all unexpected
 *   dependencies here.
 *
 */
  
/* Include guards */
#pragma once
#ifndef UTILS_H
#define UTILS_H
 
/* Forward declarations of classes that have interdependencies. */
class DRW;
class RabinPair;
class NBW;
class SafraNode;
class SafraTree;
class Conjunction;
class Literal;
 
/* Necessary to include these here for the definition of custom types.
 */
#include <boost/dynamic_bitset/dynamic_bitset.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/graph/adjacency_list.hpp>
 
/* Whether we should attempt to cache transitions of the Buchi automata.
 * If set to "true", caching procedure still depends on NBW_MAX_CACHED_SIZE
 */
#define NBW_USE_CACHE false
 
/** The maximum size an NBW can be for us to attempt to cache its states.
 * Setting this above 32 is instant suicide and setting it anywhere above 10
 * or so will be extremely prohibitive.
 */
#define NBW_MAX_CACHED_SIZE 10

/**
 * Converts an int to a string using boost. Otherwise you can accidentally 
 * append characters to strings when you're dealing with ints in the ASCII
 * char range (ecas!)
 */
#define INT_TO_STR(x) (boost::lexical_cast<std::string>(x))

/** Boundary conditions for a cellular automaton. Omega is one-way-infinite;
 * zeta is two-way-infinite.
 */
enum Boundary { OMEGA, ZETA };

/**
 * The identity cellular automaton.
 *
 */
#define IDENTITY_ECA_NUM (204)

/**
 * Defined in utils.cpp, this represents the default alphabet which labels an
 * automaton. The default value is a,b,c,d...
 */
extern std::string default_alphabet;

/* The adjacency list format which is used for Boost graphs looks like this */
typedef boost::adjacency_list< boost::vecS, 
                               boost::vecS, 
                               boost::directedS, 
                               boost::property<boost::vertex_name_t, int> > BoostGraph;

// Used to talk about a subset of the states of an automaton
typedef boost::dynamic_bitset<unsigned long> state_set_t;

/** Used to talk about a character in the alphabet of an automaton, composed
 * from one character from each of the CA tracks at the same position.
 */
typedef boost::dynamic_bitset<unsigned long> slice;

// TODO(jgershen): make better sense of all these typedefs.
/* My original intention was to use 'bitvector' to refer to a set of bits like this when
 * they're not a set of NBA states or a slice of characters, but I am beginning to think that 
 * all of these should be unified under the label of bitvector.
 */
typedef boost::dynamic_bitset<unsigned long> bitvector;


/** Defined in utils.cpp:
 *  Read a line from an input file, ignoring lines that begin with '#'.
 */
void get_next_line(std::istream &input, std::string &buffer);

#endif
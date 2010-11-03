/** @file SafraTest.hpp
 * Describes a number of functions for testing the correctness and performance
 * of Safra's construction, or for determinizing sample automata using this
 * construction.
 *
 */

#include "NBW.hpp"
#include "SafraTree.hpp"

int run_random_trials(int argc, char** argv);

int print_averages(int argc, char** argv);

void determinize(char* infile, char* outfile);

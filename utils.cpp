/**  @file utils.cpp
 *   @author Joe Gershenson
 *   
 *   Utility functions and fields that are used in multiple files or classes.
 *
 */
 
/* Include string IO + stream utilities for tokenizing input files */
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream> 
 
#include "utils.hpp"

std::string default_alphabet = "abcdefghijklmnopqrstuvwxyz";

/** Read a line from an input file, ignoring lines that begin with '#'.
 */
void get_next_line(std::istream &input, std::string &buffer){
    using namespace std;
    while( !input.eof() && input.good() ){
        getline(input, buffer);
        if(buffer[0] != '#')
            break;
    }
}
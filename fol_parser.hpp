#pragma once
#ifndef FOL_PARSER_H
#define FOL_PARSER_H

#include <string>
#include <iostream>
#include <sstream>

class InputWrapper{
    private:
        static std::stringstream* content_stream;
        static bool has_content;

    public:
        static void set_input(const std::string& input);
        static char get_char();
        static void close();
        
        static int run_parser();
};


#endif
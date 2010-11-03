#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include "NBW.hpp"
#include "logic.hpp"
#include "fol_parser.hpp"
#include "arg_parser.hpp"

namespace {

const char * invocation_name = 0;
const char * const Program_name    = "cave";
const char * const program_name    = "cave";
const char * const program_year    = "2010";

void show_help( const bool verbose ){
    std::printf( "%s - Cellular Automata Verification Environment.\n", Program_name );
    std::printf( "If you have a better acronym, let me know. At least this is a word.\n");
    std::printf( "Returns >0 if formula is valid, 0 if formula false, <0 on syntax errors.\n");
    std::printf( "See http://tenji.cdm.cs.cmu.edu/ for details.\n");
    std::printf( "\nUsage: %s [options]\n", invocation_name );
    std::printf( "\nOptions:\n" );
    std::printf( "  -h, --help                   display this help and exit\n" );
    std::printf( "  -V, --version                output version information and exit\n" );
    std::printf( "  -e, --eca=<n>                set model to use ECA n.\n" );
    std::printf( "  -f, --formula=\"<arg>\"        parse the formula instead of reading from stdin\n" );
    std::printf( "  -Z, --zeta                   work with bi-infinite cellular automata (EXPERIMENTAL)\n" );    
    std::printf( "  -v, --verbose                verbose mode\n" );
}


void show_version(){
    std::printf( "%s\n", Program_name);
    std::printf( "Copyright (C) %s Joseph Gershenson.\n", program_year );
    std::printf( "License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>\n" );
    std::printf( "This is free software: you are free to change and redistribute it.\n" );
    std::printf( "There is NO WARRANTY, to the extent permitted by law.\n" );
}


void show_error( const char * msg, const int errcode = 0, const bool help = false ){
    if( msg && msg[0] != 0 ){
        std::fprintf( stderr, "%s", msg );
        if( errcode > 0 ) 
            std::fprintf( stderr, ": %s", std::strerror( errcode ) );
        std::fprintf( stderr, "\n" );
    }
    if( help && invocation_name && invocation_name[0] != 0 )
        std::fprintf( stderr, "Try `%s --help' for more information.\n", invocation_name );
}


void internal_error( const char * msg ){
  char buf[80];
  std::snprintf( buf, sizeof( buf ), "internal error: %s.\n", msg );
  show_error( buf );
  std::exit( -3 );
}


const char * optname( const int code, const Arg_parser::Option options[] ){
    static char buf[2] = "?";
    if( code != 0 ){
        for( int i = 0; options[i].code; ++i ){
            if( code == options[i].code ){ 
                if( options[i].name ) 
                    return options[i].name; 
                else 
                    break; 
            }
        }
    }
    if( code > 0 && code < 256 ) 
        buf[0] = code; 
    else 
        buf[0] = '?';
    return buf;
}

} // end namespace


int main( int argc, char **argv)
{
    bool verbose = false;
    invocation_name = argv[0];
    
    Boundary conditions = OMEGA;

    std::string formula;
    
    
    const Arg_parser::Option options[] =
    {
        { 'V', "version",  Arg_parser::no    },
        { 'e', "eca",    Arg_parser::yes   },
        { 'f', "formula",   Arg_parser::yes },
        { 'h', "help",     Arg_parser::no    },
        { 'Z', "zeta",  Arg_parser::no    },
        { 'v', "verbose",  Arg_parser::no    },
        { 256, "orphan",   Arg_parser::no    },
        {   0, 0,          Arg_parser::no    } 
    };
    
    Arg_parser parser( argc, argv, options );
    if( parser.error().size() ){
        show_error( parser.error().c_str(), 0, true ); 
        return -1; 
    }      
    
    for( int i = 0; i < parser.arguments(); i++ ){
        const int code = parser.code(i);
        if( !code ) break;					// no more options
        switch( code ){
            case 'V': 
                show_version(); 
                return 0;
            case 'e':  // set eca
            {
                Literal::set_default_eca(atoi(parser.argument(i).c_str()));
                break;			
            }
            case 'f':  // set formula
            {
                // Add newline character to the end of the formula if neccesary
                std::string temp(parser.argument(i));
                if(temp.at(temp.length()-1) != '\n')
                    temp.push_back('\n');
                formula = temp;
                break;
            }
            case 'h': show_help( verbose ); return 0;
            case 'Z': conditions = ZETA; break;
            case 'v': verbose = true; break;
            case 256: break;				// example, do nothing
            default : internal_error( "uncaught option" );
        }
    } // end process options
    
    
    
    
    if( formula.length() > 0 ) {
        InputWrapper::set_input(formula);
    } else {
        std::cout << ">> ";
    }
    
    int parse_errcode = InputWrapper::run_parser();
    if(parse_errcode){
        fprintf(stderr, "Error parsing formula: %s\n", formula.c_str());    
        std::exit(-1);
    }
        

    std::vector<Conjunction*>* result = Conjunction::last_formula_parsed;
    for(int i = 0; i < result->size(); i++) 
        std::cout << result->at(i)->to_string() << std::endl;    



    /* Check to make sure that the user is not trying to negate a zeta-automaton.
     * This requires implementing a procedure in
     *          Cellular automata, ωω-regular sets, and sofic systems
     *              by Culik and Yu. See Discrete Applied Mathematics, 
     *              Volume 32, Issue 2, 29 July 1991, Pages 85-101 
     * Because of time constraints this is not yet a feature. 
     */
    if(conditions == ZETA){
        for(int i = 0; i < result->size(); i++) {
            for(int j = 0; j < result->at(i)->quantifiers.size(); j++) {
                if(result->at(i)->quantifiers[j].negated){
                    fprintf(stderr, "Error: negation not supported for zeta-automata.\n");
                    std::exit(-1);
                }
            }
        }
        
    }



    NBW* nbw = NBW::build_automaton(*result, conditions);
    int valid = nbw->is_empty() ? 0 : 1; // the formula is valid if no counterexamples can be found
    if( valid )
        printf("true\n");
    else
        printf("false\n");
    
    delete nbw;

    return valid;  
}

# libraries... boost etc.
boost = /usr/local/boost_1_40_0

# I will accept having to rebuild a ton of things whenever part of the spec changes.
headers = buchi_gen.hpp logic.hpp SafraTest.hpp SafraTree.hpp NBW.hpp DRW.hpp utils.hpp arg_parser.hpp

shared_objects =  NBW.o DRW.o utils.o SafraTree.o buchi_gen.o logic.o 
safra_objects = SafraTest.o 
bgen_objects = gen_test.o 
io_objects = cli.o arg_parser.o fol_parser.o

# targets are for cleanup purposes
targets = safra bgen cave

# set to -pg to enable profiling
prof_flags = 
# set to -O2 or -O3 to enable optimization
OPTIMIZE = -O3

# set compiler
CC = g++
CXX = g++

# set flags
CPPFLAGS = $(OPTIMIZE) $(prof_flags) -g -fopenmp -I$(boost) -Wno-deprecated
LDFLAGS = $(OPTIMIZE) $(prof_flags) -g -fopenmp

all: cave safra bgen

cave: $(shared_objects) $(io_objects) $(headers) utils.hpp
	g++ $(LDFLAGS) -o cave $(shared_objects) $(io_objects) -I$(boost)

safra: $(safra_objects) $(shared_objects) $(headers) utils.hpp
	g++ $(LDFLAGS) -o safra $(safra_objects) $(shared_objects) -I$(boost)

bgen: $(bgen_objects) $(shared_objects) $(headers) utils.hpp
	g++ $(LDFLAGS) -o bgen $(bgen_objects) $(shared_objects) -I$(boost)

clean:
	-rm *~ *.o $(targets)

fol_parser.o: 	fol.ypp
	bison -vt fol.ypp -o fol_parser.cpp
	g++ $(CPPFLAGS) -c -o fol_parser.o fol_parser.cpp
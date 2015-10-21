# Change this if the default compiler does not work.
#CXX = g++

# Use 64-bit integers in a 64-bit environment.
SIZE_FLAGS = -DMASSIVE_DATA_RLCSA

# Parallelism is supported by either libstdc++ Parallel Mode or MCSTL.
PARALLEL_FLAGS = #-DMULTITHREAD_SUPPORT -D_GLIBCXX_PARALLEL -fopenmp
# MCSTL_ROOT = /fs-3/d/jltsiren/suds/mcstl
# PARALLEL_FLAGS = -DMULTITHREAD_SUPPORT -I$(MCSTL_ROOT)/c++ -fopenmp

# Vectors using nibble codes instead of delta codes are faster, but they also
# take up more space.
VECTOR_FLAGS = $(PSI_FLAGS) $(LCP_FLAGS) $(SA_FLAGS)
# PSI_FLAGS = -DUSE_NIBBLE_VECTORS
# LCP_FLAGS = -DSUCCINCT_LCP_VECTOR
# SA_FLAGS = -DSUCCINCT_SA_VECTOR

RLCSA_ROOT =..
RLCSA = $(RLCSA_ROOT)/rlcsa.a

# DEBUG_FLAGS = -DDEBUG
SDSL_PREFIX = /s/chopin/l/grad/muggli/local
SDSL_STUFF =  -I/usr/include/boost/ -I$(SDSL_PREFIX)/include -L /usr/lib64 -L$(SDSL_PREFIX)/lib  -lsdsl -ldivsufsort -ldivsufsort64 

CXXFLAGS = -Wall -O3 -march=native $(DEBUG_FLAGS) $(SIZE_FLAGS) $(PARALLEL_FLAGS) $(VECTOR_FLAGS) -I$(RLCSA_ROOT) $(SDSL_STUFF) -std=c++11 

OBJS = graph.o gcsa.o parameter_handler.o
LIBRARIES = gcsa.a $(RLCSA) ~/local/lib/libsdsl.a ~/local/lib/libdivsufsort.a  ~/local/lib/libdivsufsort64.a
PROGRAMS = build_automaton clean_alignment generator build_index gcsa_test transpose  verify determinize # rlcsa_test

all: default

default: gcsa.a $(PROGRAMS)


depend:
	g++ -MM *.cpp > dependencies.mk

gcsa.a: $(OBJS)
	ar rcs gcsa.a $(OBJS)

build_automaton: build_automaton.o $(LIBRARIES)
	$(CXX) $(CXXFLAGS) -o build_automaton build_automaton.o $(LIBRARIES) 

clean_alignment: clean_alignment.o $(RLCSA)
	$(CXX) $(CXXFLAGS) -o clean_alignment clean_alignment.o $(RLCSA)

transpose: transpose.o $(RLCSA)
	$(CXX) $(CXXFLAGS) -o transpose transpose.o $(RLCSA)

verify: verify.o $(LIBRARIES)
	$(CXX) $(CXXFLAGS) -o verify verify.o $(LIBRARIES)

determinize: determinize.o $(LIBRARIES)
	$(CXX) $(CXXFLAGS) -o $@ $@.o $(LIBRARIES)

generator: generator.o
	$(CXX) $(CXXFLAGS) -o generator generator.o

build_index: build_index.o $(LIBRARIES)
	$(CXX) $(CXXFLAGS) -o build_index build_index.o $(LIBRARIES)

gcsa_test: gcsa_test.o gcsa.a $(LIBRARIES)
	$(CXX) $(CXXFLAGS) -o gcsa_test gcsa_test.o $(LIBRARIES) ../om_set1/scoringlib.o

rlcsa_test: rlcsa_test.o gcsa.a $(LIBRARIES)
	$(CXX) $(CXXFLAGS) -o rlcsa_test rlcsa_test.o $(LIBRARIES)

clean:
	rm -f *.o gcsa.a
	rm -f $(PROGRAMS)

package:
	mkdir gcsa
	cp LICENSE Makefile README dependencies.mk *.cpp *.h targz gcsa
	./targz gcsa
	rm -r gcsa/*
	rmdir gcsa

deploy:
	find . -type d -exec chmod 770 {} \;
	find . -type f -not \( -perm +111 \) -exec chmod 660 {} \;
	find . -type f -perm +111 -exec chmod 770 {} \;


include dependencies.mk

# Change this if the default compiler does not work.
CXX = g++-mp-4.8

# Use 64-bit integers in a 64-bit environment.
SIZE_FLAGS = -DMASSIVE_DATA_RLCSA

# Parallelism is supported by either libstdc++ Parallel Mode or MCSTL.
PARALLEL_FLAGS = -DMULTITHREAD_SUPPORT -D_GLIBCXX_PARALLEL -fopenmp
# MCSTL_ROOT = /fs-3/d/jltsiren/suds/mcstl
# PARALLEL_FLAGS = -DMULTITHREAD_SUPPORT -I$(MCSTL_ROOT)/c++ -fopenmp

# Vectors using nibble codes instead of delta codes are faster, but they also
# take up more space.
VECTOR_FLAGS = $(PSI_FLAGS) $(LCP_FLAGS) $(SA_FLAGS)
# PSI_FLAGS = -DUSE_NIBBLE_VECTORS
# LCP_FLAGS = -DSUCCINCT_LCP_VECTOR
# SA_FLAGS = -DSUCCINCT_SA_VECTOR

RLCSA_ROOT =../rlcsa/
RLCSA = $(RLCSA_ROOT)/rlcsa.a

GCSA_ROOT =../gcsa/
GCSA = $(GCSA_ROOT)/gcsa.a

# DEBUG_FLAGS = -DDEBUG

CXXFLAGS = -Wall -O3 $(DEBUG_FLAGS) $(SIZE_FLAGS) $(PARALLEL_FLAGS) $(VECTOR_FLAGS) -I$(RLCSA_ROOT) -I$(GCSA_ROOT)

LIBRARIES = $(GCSA) $(RLCSA)
PROGRAMS = gcsa_alignment materescue seedrescue2materescue seedrescue vcf2malignment


default: $(PROGRAMS)

gcsa_alignment: gcsa_alignment.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<  $(LIBRARIES)

materescue: materescue.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<  $(LIBRARIES)

seedrescue2materescue: seedrescue2materescue.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<  $(LIBRARIES)

seedrescue: seedrescue.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<  $(LIBRARIES)

vcf2malignment: vcf2malignment.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<  $(LIBRARIES)

clean:
	rm -f *.o
	rm -f $(PROGRAMS)

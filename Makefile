#
# Copyright 2014, Daehwan Kim <infphilo@gmail.com>
#
# This file is part of HISAT, which is copied and modified from Makefile in the Bowtie2 package.
#
# HISAT is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# HISAT is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with HISAT.  If not, see <http://www.gnu.org/licenses/>.
#
#
# Makefile for hisat-align, hisat-build, hisat-inspect
#

INC =
GCC_PREFIX = $(shell dirname `which gcc`)
GCC_SUFFIX =
CC = $(GCC_PREFIX)/gcc$(GCC_SUFFIX)
CPP = $(GCC_PREFIX)/g++$(GCC_SUFFIX)
CXX = $(CPP)
HEADERS = $(wildcard *.h)
BOWTIE_MM = 1
BOWTIE_SHARED_MEM = 0

# Detect Cygwin or MinGW
WINDOWS = 0
CYGWIN = 0
MINGW = 0
ifneq (,$(findstring CYGWIN,$(shell uname)))
	WINDOWS = 1 
	CYGWIN = 1
	# POSIX memory-mapped files not currently supported on Windows
	BOWTIE_MM = 0
	BOWTIE_SHARED_MEM = 0
else
	ifneq (,$(findstring MINGW,$(shell uname)))
		WINDOWS = 1
		MINGW = 1
		# POSIX memory-mapped files not currently supported on Windows
		BOWTIE_MM = 0
		BOWTIE_SHARED_MEM = 0
	endif
endif

MACOS = 0
ifneq (,$(findstring Darwin,$(shell uname)))
	MACOS = 1
endif

POPCNT_CAPABILITY ?= 1
ifeq (1, $(POPCNT_CAPABILITY))
    EXTRA_FLAGS += -DPOPCNT_CAPABILITY
    INC += -I third_party
endif

MM_DEF = 

ifeq (1,$(BOWTIE_MM))
	MM_DEF = -DBOWTIE_MM
endif

SHMEM_DEF = 

ifeq (1,$(BOWTIE_SHARED_MEM))
	SHMEM_DEF = -DBOWTIE_SHARED_MEM
endif

PTHREAD_PKG =
PTHREAD_LIB = 

ifeq (1,$(MINGW))
	PTHREAD_LIB = 
else
	PTHREAD_LIB = -lpthread
endif

SEARCH_LIBS = 
BUILD_LIBS = 
INSPECT_LIBS =

ifeq (1,$(MINGW))
	BUILD_LIBS = 
	INSPECT_LIBS = 
endif

USE_SRA = 0
SRA_DEF =
SRA_LIB =
SERACH_INC = 
ifeq (1,$(USE_SRA))
	SRA_DEF = -DUSE_SRA
	SRA_LIB = -lncbi-ngs-c++-static -lngs-c++-static -lncbi-vdb-static -ldl
	SEARCH_INC += -I$(NCBI_NGS_DIR)/include -I$(NCBI_VDB_DIR)/include
	SEARCH_LIBS += -L$(NCBI_NGS_DIR)/lib64 -L$(NCBI_VDB_DIR)/lib64
endif

LIBS = $(PTHREAD_LIB)

SHARED_CPPS = ccnt_lut.cpp ref_read.cpp alphabet.cpp shmem.cpp \
	edit.cpp bt2_idx.cpp \
	reference.cpp ds.cpp multikey_qsort.cpp limit.cpp \
	random_source.cpp tinythread.cpp
SEARCH_CPPS = qual.cpp pat.cpp sam.cpp \
	read_qseq.cpp aligner_seed_policy.cpp \
	aligner_seed.cpp \
	aligner_seed2.cpp \
	aligner_sw.cpp \
	aligner_sw_driver.cpp aligner_cache.cpp \
	aligner_result.cpp ref_coord.cpp mask.cpp \
	pe.cpp aln_sink.cpp dp_framer.cpp \
	scoring.cpp presets.cpp unique.cpp \
	simple_func.cpp \
	random_util.cpp \
	aligner_bt.cpp sse_util.cpp \
	aligner_swsse.cpp outq.cpp \
	aligner_swsse_loc_i16.cpp \
	aligner_swsse_ee_i16.cpp \
	aligner_swsse_loc_u8.cpp \
	aligner_swsse_ee_u8.cpp \
	aligner_driver.cpp \
	splice_site.cpp 

BUILD_CPPS = diff_sample.cpp

HISAT_CPPS_MAIN = $(SEARCH_CPPS) hisat_main.cpp
HISAT_BUILD_CPPS_MAIN = $(BUILD_CPPS) hisat_build_main.cpp

HISAT2_BUILD_CPPS_MAIN = $(BUILD_CPPS) hisat2_build_main.cpp

SEARCH_FRAGMENTS = $(wildcard search_*_phase*.c)
VERSION = $(shell cat VERSION)

# Convert BITS=?? to a -m flag
BITS=32
ifeq (x86_64,$(shell uname -m))
BITS=64
endif
# msys will always be 32 bit so look at the cpu arch instead.
ifneq (,$(findstring AMD64,$(PROCESSOR_ARCHITEW6432)))
	ifeq (1,$(MINGW))
		BITS=64
	endif
endif
BITS_FLAG =

ifeq (32,$(BITS))
	BITS_FLAG = -m32
endif

ifeq (64,$(BITS))
	BITS_FLAG = -m64
endif
SSE_FLAG=-msse2

DEBUG_FLAGS    = -O0 -g3 $(BIToS_FLAG) $(SSE_FLAG)
DEBUG_DEFS     = -DCOMPILER_OPTIONS="\"$(DEBUG_FLAGS) $(EXTRA_FLAGS)\""
RELEASE_FLAGS  = -O3 $(BITS_FLAG) $(SSE_FLAG) -funroll-loops -g3
RELEASE_DEFS   = -DCOMPILER_OPTIONS="\"$(RELEASE_FLAGS) $(EXTRA_FLAGS)\""
NOASSERT_FLAGS = -DNDEBUG
FILE_FLAGS     = -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -D_GNU_SOURCE

ifeq (1,$(USE_SRA))
	ifeq (1, $(MACOS))
		DEBUG_FLAGS += -mmacosx-version-min=10.6
		RELEASE_FLAGS += -mmacosx-version-min=10.6
	endif
endif


HISAT_BIN_LIST = hisat-build-s \
	hisat-build-l \
	hisat-align-s \
	hisat-align-l \
	hisat-inspect-s \
	hisat-inspect-l
HISAT_BIN_LIST_AUX = hisat-build-s-debug \
	hisat-build-l-debug \
	hisat-align-s-debug \
	hisat-align-l-debug \
	hisat-inspect-s-debug \
	hisat-inspect-l-debug

GENERAL_LIST = $(wildcard scripts/*.sh) \
	$(wildcard scripts/*.pl) \
	$(wildcard *.py) \
	doc/manual.inc.html \
	doc/README \
	doc/style.css \
	$(wildcard example/index/*.bt2) \
	$(wildcard example/reads/*.fq) \
	$(wildcard example/reads/*.pl) \
	example/reference/22_20-21M.fa \
	$(PTHREAD_PKG) \
	hisat \
	hisat-build \
	hisat-inspect \
	AUTHORS \
	LICENSE \
	NEWS \
	MANUAL \
	MANUAL.markdown \
	TUTORIAL \
	VERSION

ifeq (1,$(WINDOWS))
	HISAT_BIN_LIST := $(HISAT_BIN_LIST) hisat.bat hisat-build.bat hisat-inspect.bat 
endif

# This is helpful on Windows under MinGW/MSYS, where Make might go for
# the Windows FIND tool instead.
FIND=$(shell which find)

SRC_PKG_LIST = $(wildcard *.h) \
	$(wildcard *.hh) \
	$(wildcard *.c) \
	$(wildcard *.cpp) \
	doc/strip_markdown.pl \
	Makefile \
	$(GENERAL_LIST)

BIN_PKG_LIST = $(GENERAL_LIST)

.PHONY: all allall both both-debug

all: $(HISAT_BIN_LIST)

allall: $(HISAT_BIN_LIST) $(HISAT_BIN_LIST_AUX)

both: hisat-align-s hisat-align-l  hisat-build-s hisat-build-l

both-debug: hisat-align-s-debug hisat-align-l-debug hisat-build-s-debug hisat-build-l-debug

DEFS=-fno-strict-aliasing \
     -DHISAT_VERSION="\"`cat VERSION`\"" \
     -DHISAT2_VERSION="\"`cat VERSION2`\"" \
     -DBUILD_HOST="\"`hostname`\"" \
     -DBUILD_TIME="\"`date`\"" \
     -DCOMPILER_VERSION="\"`$(CXX) -v 2>&1 | tail -1`\"" \
     $(FILE_FLAGS) \
     $(PREF_DEF) \
     $(MM_DEF) \
     $(SHMEM_DEF)

#
# hisat-build targets
#

hisat-build-s: hisat_build.cpp $(SHARED_CPPS) $(HEADERS)
	$(CXX) $(RELEASE_FLAGS) $(RELEASE_DEFS) $(EXTRA_FLAGS) \
	$(DEFS) -DBOWTIE2 $(NOASSERT_FLAGS) -Wall \
	$(INC) \
	-o $@ $< \
	$(SHARED_CPPS) $(HISAT_BUILD_CPPS_MAIN) \
	$(LIBS) $(BUILD_LIBS)

hisat-build-l: hisat_build.cpp $(SHARED_CPPS) $(HEADERS)
	$(CXX) $(RELEASE_FLAGS) $(RELEASE_DEFS) $(EXTRA_FLAGS) \
	$(DEFS) -DBOWTIE2 -DBOWTIE_64BIT_INDEX $(NOASSERT_FLAGS) -Wall \
	$(INC) \
	-o $@ $< \
	$(SHARED_CPPS) $(HISAT_BUILD_CPPS_MAIN) \
	$(LIBS) $(BUILD_LIBS)

hisat-build-s-debug: hisat_build.cpp $(SHARED_CPPS) $(HEADERS)
	$(CXX) $(DEBUG_FLAGS) $(DEBUG_DEFS) $(EXTRA_FLAGS) \
	$(DEFS) -DBOWTIE2 -Wall \
	$(INC) \
	-o $@ $< \
	$(SHARED_CPPS) $(HISAT_BUILD_CPPS_MAIN) \
	$(LIBS) $(BUILD_LIBS)

hisat-build-l-debug: hisat_build.cpp $(SHARED_CPPS) $(HEADERS)
	$(CXX) $(DEBUG_FLAGS) $(DEBUG_DEFS) $(EXTRA_FLAGS) \
	$(DEFS) -DBOWTIE2 -DBOWTIE_64BIT_INDEX -Wall \
	$(INC) \
	-o $@ $< \
	$(SHARED_CPPS) $(HISAT_BUILD_CPPS_MAIN) \
	$(LIBS) $(BUILD_LIBS)

#
# hisat targets
#

hisat-align-s: hisat.cpp $(SEARCH_CPPS) $(SHARED_CPPS) $(HEADERS) $(SEARCH_FRAGMENTS)
	$(CXX) $(RELEASE_FLAGS) $(RELEASE_DEFS) $(EXTRA_FLAGS) \
	$(DEFS) $(SRA_DEF) -DBOWTIE2 $(NOASSERT_FLAGS) -Wall \
	$(INC) $(SEARCH_INC) \
	-o $@ $< \
	$(SHARED_CPPS) $(HISAT_CPPS_MAIN) \
	$(LIBS) $(SRA_LIB) $(SEARCH_LIBS)

hisat-align-l: hisat.cpp $(SEARCH_CPPS) $(SHARED_CPPS) $(HEADERS) $(SEARCH_FRAGMENTS)
	$(CXX) $(RELEASE_FLAGS) $(RELEASE_DEFS) $(EXTRA_FLAGS) \
	$(DEFS) $(SRA_DEF) -DBOWTIE2 -DBOWTIE_64BIT_INDEX $(NOASSERT_FLAGS) -Wall \
	$(INC) $(SEARCH_INC) \
	-o $@ $< \
	$(SHARED_CPPS) $(HISAT_CPPS_MAIN) \
	$(LIBS) $(SRA_LIB) $(SEARCH_LIBS)

hisat-align-s-debug: hisat.cpp $(SEARCH_CPPS) $(SHARED_CPPS) $(HEADERS) $(SEARCH_FRAGMENTS)
	$(CXX) $(DEBUG_FLAGS) \
	$(DEBUG_DEFS) $(EXTRA_FLAGS) \
	$(DEFS) $(SRA_DEF) -DBOWTIE2 -Wall \
	$(INC) $(SEARCH_INC) \
	-o $@ $< \
	$(SHARED_CPPS) $(HISAT_CPPS_MAIN) \
	$(LIBS) $(SRA_LIB) $(SEARCH_LIBS)

hisat-align-l-debug: hisat.cpp $(SEARCH_CPPS) $(SHARED_CPPS) $(HEADERS) $(SEARCH_FRAGMENTS)
	$(CXX) $(DEBUG_FLAGS) \
	$(DEBUG_DEFS) $(EXTRA_FLAGS) \
	$(DEFS) $(SRA_DEF) -DBOWTIE2 -DBOWTIE_64BIT_INDEX -Wall \
	$(INC) $(SEARCH_INC) \
	-o $@ $< \
	$(SHARED_CPPS) $(HISAT_CPPS_MAIN) \
	$(LIBS) $(SRA_LIB) $(SEARCH_LIBS)

#
# hisat-inspect targets
#

hisat-inspect-s: hisat_inspect.cpp $(HEADERS) $(SHARED_CPPS)
	$(CXX) $(RELEASE_FLAGS) \
	$(RELEASE_DEFS) $(EXTRA_FLAGS) \
	$(DEFS) -DBOWTIE2 -DHISAT_INSPECT_MAIN -Wall \
	$(INC) -I . \
	-o $@ $< \
	$(SHARED_CPPS) \
	$(LIBS) $(INSPECT_LIBS)

hisat-inspect-l: hisat_inspect.cpp $(HEADERS) $(SHARED_CPPS)
	$(CXX) $(RELEASE_FLAGS) \
	$(RELEASE_DEFS) $(EXTRA_FLAGS) \
	$(DEFS) -DBOWTIE2 -DBOWTIE_64BIT_INDEX -DHISAT_INSPECT_MAIN -Wall \
	$(INC) -I . \
	-o $@ $< \
	$(SHARED_CPPS) \
	$(LIBS) $(INSPECT_LIBS)

hisat-inspect-s-debug: hisat_inspect.cpp $(HEADERS) $(SHARED_CPPS) 
	$(CXX) $(DEBUG_FLAGS) \
	$(DEBUG_DEFS) $(EXTRA_FLAGS) \
	$(DEFS) -DBOWTIE2 -DHISAT_INSPECT_MAIN -Wall \
	$(INC) -I . \
	-o $@ $< \
	$(SHARED_CPPS) \
	$(LIBS) $(INSPECT_LIBS)

hisat-inspect-l-debug: hisat_inspect.cpp $(HEADERS) $(SHARED_CPPS) 
	$(CXX) $(DEBUG_FLAGS) \
	$(DEBUG_DEFS) $(EXTRA_FLAGS) \
	$(DEFS) -DBOWTIE2 -DBOWTIE_64BIT_INDEX -DHISAT_INSPECT_MAIN -Wall \
	$(INC) -I . \
	-o $@ $< \
	$(SHARED_CPPS) \
	$(LIBS) $(INSPECT_LIBS)


#
# hisat2-build targets
#

hisat2-build-s: hisat2_build.cpp $(SHARED_CPPS) $(HEADERS)
	$(CXX) $(RELEASE_FLAGS) $(RELEASE_DEFS) $(EXTRA_FLAGS) \
	$(DEFS) -DBOWTIE2 $(NOASSERT_FLAGS) -Wall \
	$(INC) \
	-o $@ $< \
	$(SHARED_CPPS) $(HISAT2_BUILD_CPPS_MAIN) \
	$(LIBS) $(BUILD_LIBS)

hisat2-build-l: hisat2_build.cpp $(SHARED_CPPS) $(HEADERS)
	$(CXX) $(RELEASE_FLAGS) $(RELEASE_DEFS) $(EXTRA_FLAGS) \
	$(DEFS) -DBOWTIE2 -DBOWTIE_64BIT_INDEX $(NOASSERT_FLAGS) -Wall \
	$(INC) \
	-o $@ $< \
	$(SHARED_CPPS) $(HISAT2_BUILD_CPPS_MAIN) \
	$(LIBS) $(BUILD_LIBS)

hisat2-build-s-debug: hisat2_build.cpp $(SHARED_CPPS) $(HEADERS)
	$(CXX) $(DEBUG_FLAGS) $(DEBUG_DEFS) $(EXTRA_FLAGS) \
	$(DEFS) -DBOWTIE2 -Wall \
	$(INC) \
	-o $@ $< \
	$(SHARED_CPPS) $(HISAT2_BUILD_CPPS_MAIN) \
	$(LIBS) $(BUILD_LIBS)

hisat2-build-l-debug: hisat2_build.cpp $(SHARED_CPPS) $(HEADERS)
	$(CXX) $(DEBUG_FLAGS) $(DEBUG_DEFS) $(EXTRA_FLAGS) \
	$(DEFS) -DBOWTIE2 -DBOWTIE_64BIT_INDEX -Wall \
	$(INC) \
	-o $@ $< \
	$(SHARED_CPPS) $(HISAT2_BUILD_CPPS_MAIN) \
	$(LIBS) $(BUILD_LIBS)



#
# hisat-bp targets
#

hisat-bp-bin: hisat_bp.cpp $(SEARCH_CPPS) $(SHARED_CPPS) $(HEADERS) $(SEARCH_FRAGMENTS)
	$(CXX) $(RELEASE_FLAGS) $(RELEASE_DEFS) $(EXTRA_FLAGS) \
	$(DEFS) -DBOWTIE2 $(NOASSERT_FLAGS) -Wall \
	$(INC) \
	-o $@ $< \
	$(SHARED_CPPS) $(HISAT_CPPS_MAIN) \
	$(LIBS) $(SEARCH_LIBS)

hisat-bp-bin-debug: hisat_bp.cpp $(SEARCH_CPPS) $(SHARED_CPPS) $(HEADERS) $(SEARCH_FRAGMENTS)
	$(CXX) $(DEBUG_FLAGS) \
	$(DEBUG_DEFS) $(EXTRA_FLAGS) \
	$(DEFS) -DBOWTIE2 -Wall \
	$(INC) \
	-o $@ $< \
	$(SHARED_CPPS) $(HISAT_CPPS_MAIN) \
	$(LIBS) $(SEARCH_LIBS)

#
# rlcsa, gcsa, gcsa_test
#

SIREN_INC = -Isiren/rlcsa -Isiren/gcsa
SIREN_RELEASE_FLAGS = -O3 -DNDEBUG
SIREN_DEBUG_FLAGS = -O0 -g3
SIREN_SHARED_CPPS = siren/rlcsa/rlcsa.cpp \
	siren/rlcsa/rlcsa_builder.cpp \
	siren/rlcsa/sasamples.cpp \
	siren/rlcsa/alphabet.cpp \
	siren/rlcsa/lcpsamples.cpp \
	siren/rlcsa/sampler.cpp \
	siren/rlcsa/suffixarray.cpp \
	siren/rlcsa/adaptive_samples.cpp \
	siren/rlcsa/docarray.cpp \
	siren/rlcsa/bits/array.cpp \
	siren/rlcsa/bits/bitbuffer.cpp \
	siren/rlcsa/bits/multiarray.cpp \
	siren/rlcsa/bits/bitvector.cpp \
	siren/rlcsa/bits/deltavector.cpp \
	siren/rlcsa/bits/rlevector.cpp \
	siren/rlcsa/bits/nibblevector.cpp \
	siren/rlcsa/bits/succinctvector.cpp \
	siren/rlcsa/misc/parameters.cpp \
	siren/rlcsa/misc/utils.cpp \
	siren/gcsa/graph.cpp \
	siren/gcsa/gcsa.cpp \
	siren/gcsa/parameter_handler.cpp

clean_alignment: siren/gcsa/clean_alignment.cpp $(SIREN_SHARED_CPPS)
	$(CXX) $(SIREN_RELEASE_FLAGS) -Wall -DMASSIVE_DATA_RLCSA \
	$(SIREN_INC) \
	-o $@ $< \
	$(SIREN_SHARED_CPPS)

build_automaton: siren/gcsa/build_automaton.cpp $(SIREN_SHARED_CPPS)
	$(CXX) $(SIREN_RELEASE_FLAGS) -Wall -DMASSIVE_DATA_RLCSA \
	$(SIREN_INC) \
	-o $@ $< \
	$(SIREN_SHARED_CPPS)

build_automaton-debug: siren/gcsa/build_automaton.cpp $(SIREN_SHARED_CPPS)
	$(CXX) $(SIREN_DEBUG_FLAGS) -Wall -DMASSIVE_DATA_RLCSA \
	$(SIREN_INC) \
	-o $@ $< \
	$(SIREN_SHARED_CPPS)

determinize: siren/gcsa/determinize.cpp $(SIREN_SHARED_CPPS)
	$(CXX) $(SIREN_RELEASE_FLAGS) -Wall -DMASSIVE_DATA_RLCSA \
	$(SIREN_INC) \
	-o $@ $< \
	$(SIREN_SHARED_CPPS)

determinize-debug: siren/gcsa/determinize.cpp $(SIREN_SHARED_CPPS)
	$(CXX) $(SIREN_DEBUG_FLAGS) -Wall -DMASSIVE_DATA_RLCSA \
	$(SIREN_INC) \
	-o $@ $< \
	$(SIREN_SHARED_CPPS)

build_index: siren/gcsa/build_index.cpp $(SIREN_SHARED_CPPS)
	$(CXX) $(SIREN_RELEASE_FLAGS) -Wall -DMASSIVE_DATA_RLCSA \
	$(SIREN_INC) \
	-o $@ $< \
	$(SIREN_SHARED_CPPS)

build_index-debug: siren/gcsa/build_index.cpp $(SIREN_SHARED_CPPS)
	$(CXX) $(SIREN_DEBUG_FLAGS) -Wall -DMASSIVE_DATA_RLCSA \
	$(SIREN_INC) \
	-o $@ $< \
	$(SIREN_SHARED_CPPS)

gcsa_test: siren/gcsa/gcsa_test.cpp $(SIREN_SHARED_CPPS)
	$(CXX) $(SIREN_RELEASE_FLAGS) -Wall -DMASSIVE_DATA_RLCSA \
	$(SIREN_INC) \
	-o $@ $< \
	$(SIREN_SHARED_CPPS)

gcsa_test-debug: siren/gcsa/gcsa_test.cpp $(SIREN_SHARED_CPPS)
	$(CXX) $(SIREN_DEBUG_FLAGS) -Wall -DMASSIVE_DATA_RLCSA \
	$(SIREN_INC) \
	-o $@ $< \
	$(SIREN_SHARED_CPPS)

gcsa_alignment: siren/gcsa_util/gcsa_alignment.cpp $(SIREN_SHARED_CPPS)
	$(CXX) $(SIREN_RELEASE_FLAGS) -Wall -DMASSIVE_DATA_RLCSA \
	$(SIREN_INC) \
	-o $@ $< \
	$(SIREN_SHARED_CPPS)

gcsa_alignment-debug: siren/gcsa_util/gcsa_alignment.cpp $(SIREN_SHARED_CPPS)
	$(CXX) $(SIREN_DEBUG_FLAGS) -Wall -DMASSIVE_DATA_RLCSA \
	$(SIREN_INC) \
	-o $@ $< \
	$(SIREN_SHARED_CPPS)





hisat: ;

hisat.bat:
	echo "@echo off" > hisat.bat
	echo "perl %~dp0/hisat %*" >> hisat.bat

hisat-build.bat:
	echo "@echo off" > hisat-build.bat
	echo "python %~dp0/hisat-build %*" >> hisat-build.bat

hisat-inspect.bat:
	echo "@echo off" > hisat-inspect.bat
	echo "python %~dp0/hisat-inspect %*" >> hisat-inspect.bat


.PHONY: hisat-src
hisat-src: $(SRC_PKG_LIST)
	chmod a+x scripts/*.sh scripts/*.pl
	mkdir .src.tmp
	mkdir .src.tmp/hisat-$(VERSION)
	zip tmp.zip $(SRC_PKG_LIST)
	mv tmp.zip .src.tmp/hisat-$(VERSION)
	cd .src.tmp/hisat-$(VERSION) ; unzip tmp.zip ; rm -f tmp.zip
	cd .src.tmp ; zip -r hisat-$(VERSION)-source.zip hisat-$(VERSION)
	cp .src.tmp/hisat-$(VERSION)-source.zip .
	rm -rf .src.tmp

.PHONY: hisat-bin
hisat-bin: $(BIN_PKG_LIST) $(HISAT_BIN_LIST) $(HISAT_BIN_LIST_AUX)
	chmod a+x scripts/*.sh scripts/*.pl
	rm -rf .bin.tmp
	mkdir .bin.tmp
	mkdir .bin.tmp/hisat-$(VERSION)
	if [ -f hisat.exe ] ; then \
		zip tmp.zip $(BIN_PKG_LIST) $(addsuffix .exe,$(HISAT_BIN_LIST) $(HISAT_BIN_LIST_AUX)) ; \
	else \
		zip tmp.zip $(BIN_PKG_LIST) $(HISAT_BIN_LIST) $(HISAT_BIN_LIST_AUX) ; \
	fi
	mv tmp.zip .bin.tmp/hisat-$(VERSION)
	cd .bin.tmp/hisat-$(VERSION) ; unzip tmp.zip ; rm -f tmp.zip
	cd .bin.tmp ; zip -r hisat-$(VERSION)-$(BITS).zip hisat-$(VERSION)
	cp .bin.tmp/hisat-$(VERSION)-$(BITS).zip .
	rm -rf .bin.tmp

.PHONY: doc
doc: doc/manual.inc.html MANUAL

doc/manual.inc.html: MANUAL.markdown
	pandoc -T "HISAT Manual" -o $@ \
	 --from markdown --to HTML --toc $^
	perl -i -ne \
	 '$$w=0 if m|^</body>|;print if $$w;$$w=1 if m|^<body>|;' $@

MANUAL: MANUAL.markdown
	perl doc/strip_markdown.pl < $^ > $@

.PHONY: clean
clean:
	rm -f $(HISAT_BIN_LIST) $(HISAT_BIN_LIST_AUX) \
	$(addsuffix .exe,$(HISAT_BIN_LIST) $(HISAT_BIN_LIST_AUX)) \
	hisat-src.zip hisat-bin.zip
	rm -f core.* .tmp.head
	rm -rf *.dSYM

.PHONY: push-doc
push-doc: doc/manual.inc.html
	scp doc/*.*html igm1:/data1/igm3/www/ccb.jhu.edu/html/software/hisat/
	
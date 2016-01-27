#
# Copyright 2015, Daehwan Kim <infphilo@gmail.com>
#
# This file is part of HISAT2.
#
# HISAT 2 is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# HISAT 2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with HISAT.  If not, see <http://www.gnu.org/licenses/>.
#
#
# Makefile for hisat2-align, hisat2-build, hisat2-inspect
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

EXTRA_FLAGS += -DPOPCNT_CAPABILITY
INC += -I third_party

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
	edit.cpp gfm.cpp \
	reference.cpp ds.cpp multikey_qsort.cpp limit.cpp \
	random_source.cpp tinythread.cpp
SEARCH_CPPS = qual.cpp pat.cpp \
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

HISAT2_CPPS_MAIN = $(SEARCH_CPPS) hisat2_main.cpp
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


HISAT2_BIN_LIST = hisat2-build-s \
	hisat2-build-l \
	hisat2-align-s \
	hisat2-align-l \
	hisat2-inspect-s \
	hisat2-inspect-l
HISAT2_BIN_LIST_AUX = hisat2-build-s-debug \
	hisat2-build-l-debug \
	hisat2-align-s-debug \
	hisat2-align-l-debug \
	hisat2-inspect-s-debug \
	hisat2-inspect-l-debug

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
	hisat2 \
	hisat2-build \
	hisat2-inspect \
	AUTHORS \
	LICENSE \
	NEWS \
	MANUAL \
	MANUAL.markdown \
	TUTORIAL \
	VERSION

ifeq (1,$(WINDOWS))
	HISAT2_BIN_LIST := $(HISAT2_BIN_LIST) hisat2.bat hisat2-build.bat hisat2-inspect.bat 
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

all: $(HISAT2_BIN_LIST)

allall: $(HISAT2_BIN_LIST) $(HISAT2_BIN_LIST_AUX)

both: hisat2-align-s hisat2-align-l hisat2-build-s hisat2-build-l

both-debug: hisat2-align-s-debug hisat2-align-l-debug hisat2-build-s-debug hisat2-build-l-debug

DEFS=-fno-strict-aliasing \
     -DHISAT2_VERSION="\"`cat VERSION`\"" \
     -DBUILD_HOST="\"`hostname`\"" \
     -DBUILD_TIME="\"`date`\"" \
     -DCOMPILER_VERSION="\"`$(CXX) -v 2>&1 | tail -1`\"" \
     $(FILE_FLAGS) \
     $(PREF_DEF) \
     $(MM_DEF) \
     $(SHMEM_DEF)

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
# hisat2-build targets
#

hisat2-build-s: hisat2_build.cpp $(SHARED_CPPS) $(HEADERS)
	$(CXX) $(RELEASE_FLAGS) $(RELEASE_DEFS) $(EXTRA_FLAGS) \
	$(DEFS) -DBOWTIE2 $(NOASSERT_FLAGS) -Wall -DMASSIVE_DATA_RLCSA \
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
	$(DEFS) -DBOWTIE2 -Wall -DMASSIVE_DATA_RLCSA \
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
# hisat2 targets
#

hisat2-align-s: hisat2.cpp $(SEARCH_CPPS) $(SHARED_CPPS) $(HEADERS) $(SEARCH_FRAGMENTS)
	$(CXX) $(RELEASE_FLAGS) $(RELEASE_DEFS) $(EXTRA_FLAGS) \
	$(DEFS) $(SRA_DEF) -DBOWTIE2 $(NOASSERT_FLAGS) -Wall \
	$(INC) $(SEARCH_INC) \
	-o $@ $< \
	$(SHARED_CPPS) $(HISAT2_CPPS_MAIN) \
	$(LIBS) $(SRA_LIB) $(SEARCH_LIBS)

hisat2-align-l: hisat2.cpp $(SEARCH_CPPS) $(SHARED_CPPS) $(HEADERS) $(SEARCH_FRAGMENTS)
	$(CXX) $(RELEASE_FLAGS) $(RELEASE_DEFS) $(EXTRA_FLAGS) \
	$(DEFS) $(SRA_DEF) -DBOWTIE2 -DBOWTIE_64BIT_INDEX $(NOASSERT_FLAGS) -Wall \
	$(INC) $(SEARCH_INC) \
	-o $@ $< \
	$(SHARED_CPPS) $(HISAT2_CPPS_MAIN) \
	$(LIBS) $(SRA_LIB) $(SEARCH_LIBS)

hisat2-align-s-debug: hisat2.cpp $(SEARCH_CPPS) $(SHARED_CPPS) $(HEADERS) $(SEARCH_FRAGMENTS)
	$(CXX) $(DEBUG_FLAGS) \
	$(DEBUG_DEFS) $(EXTRA_FLAGS) \
	$(DEFS) $(SRA_DEF) -DBOWTIE2 -Wall \
	$(INC) $(SEARCH_INC) \
	-o $@ $< \
	$(SHARED_CPPS) $(HISAT2_CPPS_MAIN) \
	$(LIBS) $(SRA_LIB) $(SEARCH_LIBS)

hisat2-align-l-debug: hisat2.cpp $(SEARCH_CPPS) $(SHARED_CPPS) $(HEADERS) $(SEARCH_FRAGMENTS)
	$(CXX) $(DEBUG_FLAGS) \
	$(DEBUG_DEFS) $(EXTRA_FLAGS) \
	$(DEFS) $(SRA_DEF) -DBOWTIE2 -DBOWTIE_64BIT_INDEX -Wall \
	$(INC) $(SEARCH_INC) \
	-o $@ $< \
	$(SHARED_CPPS) $(HISAT2_CPPS_MAIN) \
	$(LIBS) $(SRA_LIB) $(SEARCH_LIBS)

#
# hisat2-inspect targets
#

hisat2-inspect-s: hisat2_inspect.cpp $(HEADERS) $(SHARED_CPPS)
	$(CXX) $(RELEASE_FLAGS) \
	$(RELEASE_DEFS) $(EXTRA_FLAGS) \
	$(DEFS) -DBOWTIE2 -DHISAT2_INSPECT_MAIN -Wall \
	$(INC) -I . \
	-o $@ $< \
	$(SHARED_CPPS) \
	$(LIBS) $(INSPECT_LIBS)

hisat2-inspect-l: hisat2_inspect.cpp $(HEADERS) $(SHARED_CPPS)
	$(CXX) $(RELEASE_FLAGS) \
	$(RELEASE_DEFS) $(EXTRA_FLAGS) \
	$(DEFS) -DBOWTIE2 -DBOWTIE_64BIT_INDEX -DHISAT2_INSPECT_MAIN -Wall \
	$(INC) -I . \
	-o $@ $< \
	$(SHARED_CPPS) \
	$(LIBS) $(INSPECT_LIBS)

hisat2-inspect-s-debug: hisat2_inspect.cpp $(HEADERS) $(SHARED_CPPS) 
	$(CXX) $(DEBUG_FLAGS) \
	$(DEBUG_DEFS) $(EXTRA_FLAGS) \
	$(DEFS) -DBOWTIE2 -DHISAT2_INSPECT_MAIN -Wall \
	$(INC) -I . \
	-o $@ $< \
	$(SHARED_CPPS) \
	$(LIBS) $(INSPECT_LIBS)

hisat2-inspect-l-debug: hisat2_inspect.cpp $(HEADERS) $(SHARED_CPPS) 
	$(CXX) $(DEBUG_FLAGS) \
	$(DEBUG_DEFS) $(EXTRA_FLAGS) \
	$(DEFS) -DBOWTIE2 -DBOWTIE_64BIT_INDEX -DHISAT2_INSPECT_MAIN -Wall \
	$(INC) -I . \
	-o $@ $< \
	$(SHARED_CPPS) \
	$(LIBS) $(INSPECT_LIBS)



hisat2: ;

hisat2.bat:
	echo "@echo off" > hisat2.bat
	echo "perl %~dp0/hisat2 %*" >> hisat2.bat

hisat2-build.bat:
	echo "@echo off" > hisat2-build.bat
	echo "python %~dp0/hisat2-build %*" >> hisat2-build.bat

hisat2-inspect.bat:
	echo "@echo off" > hisat2-inspect.bat
	echo "python %~dp0/hisat2-inspect %*" >> hisat2-inspect.bat


.PHONY: hisat2-src
hisat2-src: $(SRC_PKG_LIST)
	chmod a+x scripts/*.sh scripts/*.pl
	mkdir .src.tmp
	mkdir .src.tmp/hisat2-$(VERSION)
	zip tmp.zip $(SRC_PKG_LIST)
	mv tmp.zip .src.tmp/hisat2-$(VERSION)
	cd .src.tmp/hisat2-$(VERSION) ; unzip tmp.zip ; rm -f tmp.zip
	cd .src.tmp ; zip -r hisat2-$(VERSION)-source.zip hisat2-$(VERSION)
	cp .src.tmp/hisat2-$(VERSION)-source.zip .
	rm -rf .src.tmp

.PHONY: hisat2-bin
hisat2-bin: $(BIN_PKG_LIST) $(HISAT2_BIN_LIST) $(HISAT2_BIN_LIST_AUX)
	chmod a+x scripts/*.sh scripts/*.pl
	rm -rf .bin.tmp
	mkdir .bin.tmp
	mkdir .bin.tmp/hisat2-$(VERSION)
	if [ -f hisat2.exe ] ; then \
		zip tmp.zip $(BIN_PKG_LIST) $(addsuffix .exe,$(HISAT2_BIN_LIST) $(HISAT2_BIN_LIST_AUX)) ; \
	else \
		zip tmp.zip $(BIN_PKG_LIST) $(HISAT2_BIN_LIST) $(HISAT2_BIN_LIST_AUX) ; \
	fi
	mv tmp.zip .bin.tmp/hisat2-$(VERSION)
	cd .bin.tmp/hisat2-$(VERSION) ; unzip tmp.zip ; rm -f tmp.zip
	cd .bin.tmp ; zip -r hisat2-$(VERSION)-$(BITS).zip hisat2-$(VERSION)
	cp .bin.tmp/hisat2-$(VERSION)-$(BITS).zip .
	rm -rf .bin.tmp

.PHONY: doc
doc: doc/manual.inc.html MANUAL

doc/manual.inc.html: MANUAL.markdown
	pandoc -T "HISAT2 Manual" -o $@ \
	 --from markdown --to HTML --toc $^
	perl -i -ne \
	 '$$w=0 if m|^</body>|;print if $$w;$$w=1 if m|^<body>|;' $@

MANUAL: MANUAL.markdown
	perl doc/strip_markdown.pl < $^ > $@

.PHONY: clean
clean:
	rm -f $(HISAT2_BIN_LIST) $(HISAT2_BIN_LIST_AUX) \
	$(addsuffix .exe,$(HISAT2_BIN_LIST) $(HISAT2_BIN_LIST_AUX)) \
	hisat2-src.zip hisat2-bin.zip
	rm -f core.* .tmp.head
	rm -rf *.dSYM

.PHONY: push-doc
push-doc: doc/manual.inc.html
	scp doc/*.*html igm1:/data1/igm3/www/ccb.jhu.edu/html/software/hisat2/

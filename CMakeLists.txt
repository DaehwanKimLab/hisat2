cmake_minimum_required(VERSION 3.10)
project(hisat2)

# escape preprocessor definition 
cmake_policy(SET CMP0005 NEW)

if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()

option(BOWTIE_MM "enable memory mapping" ON)
option(BOWTIE_SHARED_MEM "enable shared memory mapping" OFF)
option(USE_SRA "enable NCBI SRA Toolkit" OFF)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_VERBOSE_MAKEFILE true)
set(CMAKE_EXPORT_COMPILE_COMMANDS on)
set(NCBI_NGS_DIR ${NCBI_NGS_DIR})

# Common Flag
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64 -g3 -Wall -msse2")

# Debug/Release Flag
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -funroll-loops")

set(INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/bin" CACHE PATH "Installation directory for executables")

if(MINGW OR CYGWIN)
    option(BOWTIE_MM "Memory mapped files not supported on Windows" OFF)
    option(BOWTIE_SHARED_MEM "Shared memory not supported on Windows" OFF)
endif()
 
if (USE_SRA AND APPLE)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mmacosx-version-min=10.6")
endif()

set(HISAT2_BIN_LIST
	hisat2-build-s
	hisat2-build-l
	hisat2-align-s
	hisat2-align-l
	hisat2-inspect-s
	hisat2-inspect-l
	hisat2-repeat
    )

set(HISAT2_LIB_LIST
	hisat2lib
	)


# Source Codes
set(SHARED_CPPS
	alphabet.cpp
	ccnt_lut.cpp
	ds.cpp
	edit.cpp
	gfm.cpp
	limit.cpp
	multikey_qsort.cpp
	random_source.cpp
	reference.cpp
	ref_read.cpp
	shmem.cpp
	tinythread.cpp
	)

set(SEARCH_CPPS
	aligner_bt.cpp
	aligner_cache.cpp
	aligner_driver.cpp
	aligner_result.cpp
	aligner_seed.cpp
	aligner_seed2.cpp
	aligner_seed_policy.cpp
	aligner_sw.cpp
	aligner_sw_driver.cpp
	aligner_swsse.cpp
	aligner_swsse_ee_i16.cpp
	aligner_swsse_ee_u8.cpp
	aligner_swsse_loc_i16.cpp
	aligner_swsse_loc_u8.cpp
	aln_sink.cpp
	dp_framer.cpp
	outq.cpp
	pat.cpp
	pe.cpp
	presets.cpp
	qual.cpp
	random_util.cpp
	read_qseq.cpp
	ref_coord.cpp mask.cpp
	scoring.cpp
	simple_func.cpp
	splice_site.cpp
	sse_util.cpp
	unique.cpp
	hisat2_main.cpp
	hisat2.cpp
	)

set(BUILD_CPPS
	diff_sample.cpp
	hisat2_build_main.cpp
	hisat2_build.cpp
	)

set(REPEAT_CPPS
	aligner_bt.cpp
	aligner_result.cpp
	aligner_sw_driver.cpp
	aligner_sw.cpp
	aligner_swsse_ee_i16.cpp
	aligner_swsse_ee_u8.cpp
	aligner_swsse_loc_i16.cpp
	aligner_swsse_loc_u8.cpp
	aligner_swsse.cpp
	bit_packed_array.cpp
	bit_packed_array.h
	dp_framer.cpp
	mask.cpp
	qual.cpp
	repeat_builder.cpp
	scoring.cpp
	simple_func.cpp
	diff_sample.cpp
	hisat2_repeat_main.cpp
	hisat2_repeat.cpp
	)

set(INSPECT_CPPS
	hisat2_inspect.cpp
	)


set(HT2LIB_CPPS
	hisat2lib/ht2_init.cpp
	hisat2lib/ht2_repeat.cpp
	hisat2lib/ht2_index.cpp
	)

set(threeN_CPPS
		alignment_3N.cpp position_3N.cpp)

string(TIMESTAMP BUILD_DATE)
file(STRINGS VERSION HISAT2_VERSION)
execute_process(COMMAND ${CMAKE_CXX_COMPILER} -v OUTPUT_VARIABLE COMPILER_VERSION)

add_definitions(
	-DBOWTIE2
    -DBUILD_HOST="${CMAKE_HOST_SYSTEM}"
	-DBUILD_TIME="${BUILD_DATE}"
    -DCOMPILER_VERSION="${CMAKE_CXX_COMPILER} ${CMAKE_CXX_COMPILER_VERSION}"
    -DHISAT2_VERSION="${HISAT2_VERSION}"
	-D_LARGEFILE_SOURCE
	-D_FILE_OFFSET_BITS=64
	-D_GNU_SOURCE
	-DPOPCNT_CAPABILITY
	)


include_directories(third_party)

link_libraries(
	pthread
	)


if(BOWTIE_MM)
	add_definitions(-DBOWTIE_MM)
endif()

if(BOWTIE_SHARED_MEM)
	add_definitions(-DBOWTIE_SHARED_MEM)
endif()

include_directories(${PROJECT_SOURCE_DIR})
add_definitions(-DCOMPILER_OPTIONS="${CMAKE_CXX_FLAGS}")

add_library(hisat2lib_static STATIC ${HT2LIB_CPPS} ${SHARED_CPPS})
add_library(hisat2lib_shared SHARED ${HT2LIB_CPPS} ${SHARED_CPPS})
add_custom_target(ht2lib)
add_dependencies(ht2lib hisat2lib_static hisat2lib_shared)

set_target_properties(hisat2lib_static PROPERTIES OUTPUT_NAME hisat2lib)
set_target_properties(hisat2lib_shared PROPERTIES OUTPUT_NAME hisat2lib)





add_executable(hisat2-align-s ${SEARCH_CPPS} ${SHARED_CPPS} ${HT2LIB_CPPS} ${threeN_CPPS})
add_executable(hisat2-align-l ${SEARCH_CPPS} ${SHARED_CPPS} ${HT2LIB_CPPS} ${threeN_CPPS})
add_executable(hisat2-build-s ${BUILD_CPPS} ${SHARED_CPPS})
add_executable(hisat2-build-l ${BUILD_CPPS} ${SHARED_CPPS})
add_executable(hisat2-inspect-s ${INSPECT_CPPS} ${SHARED_CPPS})
add_executable(hisat2-inspect-l ${INSPECT_CPPS} ${SHARED_CPPS})
add_executable(hisat2-repeat ${REPEAT_CPPS} ${SHARED_CPPS})

set_target_properties(${HISAT2_BIN_LIST} PROPERTIES DEBUG_POSTFIX "-debug")
set_target_properties(hisat2-align-l hisat2-build-l hisat2-inspect-l hisat2-repeat 
	PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS} -DBOWTIE_64BIT_INDEX")

target_compile_options(hisat2-inspect-s PUBLIC "-DHISAT2_INSPECT_MAIN")
target_compile_options(hisat2-inspect-l PUBLIC "-DHISAT2_INSPECT_MAIN")

target_compile_options(hisat2-build-s PUBLIC "-DMASSIVE_DATA_RLCSA")

target_compile_options(hisat2lib_static PUBLIC "-fPIC" "-DHISAT2_BUILD_LIB")
target_compile_options(hisat2lib_shared PUBLIC "-fPIC" "-DHISAT2_BUILD_LIB")

target_include_directories(hisat2lib_static PUBLIC hisat2lib)
target_include_directories(hisat2lib_shared PUBLIC hisat2lib)

if(USE_SRA)
	if(NOT NCBI_NGS_DIR)
		message(FATAL_ERROR "You should set NCBI_NGS_DIR path. cmake3 .... -DNCBI_NGS_DIR=/path/to/ngs")
	endif()
	#set(NCBI_NGS_DIR /home/parkch/.local/ngs)

	target_compile_options(hisat2-align-s PUBLIC "-DUSE_SRA")
	target_compile_options(hisat2-align-l PUBLIC "-DUSE_SRA")
	target_include_directories(hisat2-align-s PUBLIC "${NCBI_NGS_DIR}/include")
	target_include_directories(hisat2-align-l PUBLIC "${NCBI_NGS_DIR}/include")

	set(NCBI_LIBS 
		${NCBI_NGS_DIR}/lib64/libncbi-ngs-c++-static.a
		${NCBI_NGS_DIR}/lib64/libngs-c++-static.a
		${NCBI_NGS_DIR}/lib64/libncbi-vdb-static.a
		${CMAKE_DL_LIBS}
		)

	target_link_libraries(hisat2-align-s ${NCBI_LIBS})
	target_link_libraries(hisat2-align-l ${NCBI_LIBS})
endif()

add_custom_target(hisat2-all)
add_dependencies(hisat2-all
		hisat2-align-s
		hisat2-align-l
		hisat2-build-s
		hisat2-build-l
		hisat2-inspect-s
		hisat2-inspect-l
		hisat2-repeat)


#
# Examples
#
#target_compile_options(hisat2-inspect-s PUBLIC $<$<CONFIG:DEBUG>:-DDEBUGONLY>)
#target_compile_options(hisat2-inspect-s PUBLIC $<$<CONFIG:RELEASE>:-DRELEASEONLY>)

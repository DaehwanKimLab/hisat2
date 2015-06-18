/*
 * Copyright 2011, Ben Langmead <langmea@cs.jhu.edu>
 *
 * This file is part of Bowtie 2.
 *
 * Bowtie 2 is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Bowtie 2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Bowtie 2.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef WORD_IO_H_
#define WORD_IO_H_

#include <stdint.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include "assert_helpers.h"
#include "endian_swap.h"

/**
 * Write a 32-bit unsigned to an output stream being careful to
 * re-endianize if caller-requested endianness differs from current
 * host.
 */
static inline void writeU32(std::ostream& out, uint32_t x, bool toBigEndian) {
	uint32_t y = endianizeU32(x, toBigEndian);
	out.write((const char*)&y, 4);
}

/**
 * Write a 32-bit unsigned to an output stream using the native
 * endianness.
 */
static inline void writeU32(std::ostream& out, uint32_t x) {
	out.write((const char*)&x, 4);
}

/**
 * Write a 32-bit signed int to an output stream being careful to
 * re-endianize if caller-requested endianness differs from current
 * host.
 */
static inline void writeI32(std::ostream& out, int32_t x, bool toBigEndian) {
	int32_t y = endianizeI32(x, toBigEndian);
	out.write((const char*)&y, 4);
}

/**
 * Write a 32-bit unsigned to an output stream using the native
 * endianness.
 */
static inline void writeI32(std::ostream& out, int32_t x) {
	out.write((const char*)&x, 4);
}

/**
 * Write a 16-bit unsigned to an output stream being careful to
 * re-endianize if caller-requested endianness differs from current
 * host.
 */
static inline void writeU16(std::ostream& out, uint16_t x, bool toBigEndian) {
	uint16_t y = endianizeU16(x, toBigEndian);
	out.write((const char*)&y, 2);
}

/**
 * Write a 16-bit unsigned to an output stream using the native
 * endianness.
 */
static inline void writeU16(std::ostream& out, uint16_t x) {
	out.write((const char*)&x, 2);
}

/**
 * Write a 16-bit signed int to an output stream being careful to
 * re-endianize if caller-requested endianness differs from current
 * host.
 */
static inline void writeI16(std::ostream& out, int16_t x, bool toBigEndian) {
	int16_t y = endianizeI16(x, toBigEndian);
	out.write((const char*)&y, 2);
}

/**
 * Write a 16-bit unsigned to an output stream using the native
 * endianness.
 */
static inline void writeI16(std::ostream& out, int16_t x) {
	out.write((const char*)&x, 2);
}

/**
 * Read a 32-bit unsigned from an input stream, inverting endianness
 * if necessary.
 */
static inline uint32_t readU32(std::istream& in, bool swap) {
	uint32_t x;
	in.read((char *)&x, 4);
	assert_eq(4, in.gcount());
	if(swap) {
		return endianSwapU32(x);
	} else {
		return x;
	}
}

/**
 * Read a 32-bit unsigned from a file descriptor, optionally inverting
 * endianness.
 */
#ifdef BOWTIE_MM
static inline uint32_t readU32(int in, bool swap) {
	uint32_t x;
	if(read(in, (void *)&x, 4) != 4) {
		assert(false);
	}
	if(swap) {
		return endianSwapU32(x);
	} else {
		return x;
	}
}
#endif

/**
 * Read a 32-bit unsigned from a FILE*, optionally inverting
 * endianness.
 */
static inline uint32_t readU32(FILE* in, bool swap) {
	uint32_t x;
	if(fread((void *)&x, 1, 4, in) != 4) {
		assert(false);
	}
	if(swap) {
		return endianSwapU32(x);
	} else {
		return x;
	}
}


/**
 * Read a 32-bit signed from an input stream, inverting endianness
 * if necessary.
 */
static inline int32_t readI32(std::istream& in, bool swap) {
	int32_t x;
	in.read((char *)&x, 4);
	assert_eq(4, in.gcount());
	if(swap) {
		return endianSwapI32(x);
	} else {
		return x;
	}
}

/**
 * Read a 32-bit unsigned from a file descriptor, optionally inverting
 * endianness.
 */
#ifdef BOWTIE_MM
static inline uint32_t readI32(int in, bool swap) {
	int32_t x;
	if(read(in, (void *)&x, 4) != 4) {
		assert(false);
	}
	if(swap) {
		return endianSwapI32(x);
	} else {
		return x;
	}
}
#endif

/**
 * Read a 32-bit unsigned from a FILE*, optionally inverting
 * endianness.
 */
static inline uint32_t readI32(FILE* in, bool swap) {
	int32_t x;
	if(fread((void *)&x, 1, 4, in) != 4) {
		assert(false);
	}
	if(swap) {
		return endianSwapI32(x);
	} else {
		return x;
	}
}


/**
 * Read a 16-bit unsigned from an input stream, inverting endianness
 * if necessary.
 */
static inline uint16_t readU16(std::istream& in, bool swap) {
	uint16_t x;
	in.read((char *)&x, 2);
	assert_eq(2, in.gcount());
	if(swap) {
		return endianSwapU16(x);
	} else {
		return x;
	}
}

/**
 * Read a 16-bit unsigned from a file descriptor, optionally inverting
 * endianness.
 */
#ifdef BOWTIE_MM
static inline uint16_t readU16(int in, bool swap) {
	uint16_t x;
	if(read(in, (void *)&x, 2) != 2) {
		assert(false);
	}
	if(swap) {
		return endianSwapU16(x);
	} else {
		return x;
	}
}
#endif

/**
 * Read a 16-bit unsigned from a FILE*, optionally inverting
 * endianness.
 */
static inline uint16_t readU16(FILE* in, bool swap) {
	uint16_t x;
	if(fread((void *)&x, 1, 2, in) != 2) {
		assert(false);
	}
	if(swap) {
		return endianSwapU32(x);
	} else {
		return x;
	}
}


/**
 * Read a 16-bit signed from an input stream, inverting endianness
 * if necessary.
 */
static inline int32_t readI16(std::istream& in, bool swap) {
	int16_t x;
	in.read((char *)&x, 2);
	assert_eq(2, in.gcount());
	if(swap) {
		return endianSwapI16(x);
	} else {
		return x;
	}
}

/**
 * Read a 16-bit unsigned from a file descriptor, optionally inverting
 * endianness.
 */
#ifdef BOWTIE_MM
static inline uint16_t readI16(int in, bool swap) {
	int16_t x;
	if(read(in, (void *)&x, 2) != 2) {
		assert(false);
	}
	if(swap) {
		return endianSwapI16(x);
	} else {
		return x;
	}
}
#endif

/**
 * Read a 16-bit unsigned from a FILE*, optionally inverting
 * endianness.
 */
static inline uint16_t readI16(FILE* in, bool swap) {
	int16_t x;
	if(fread((void *)&x, 1, 2, in) != 2) {
		assert(false);
	}
	if(swap) {
		return endianSwapI16(x);
	} else {
		return x;
	}
}

template <typename index_t>
void writeIndex(std::ostream& out, index_t x, bool toBigEndian) {
	index_t y = endianizeIndex(x, toBigEndian);
	out.write((const char*)&y, sizeof(index_t));
}

/**
 * Read a unsigned from an input stream, inverting endianness
 * if necessary.
 */
template <typename index_t>
static inline index_t readIndex(std::istream& in, bool swap) {
	index_t x;
	in.read((char *)&x, sizeof(index_t));
	assert_eq(sizeof(index_t), in.gcount());
	if(swap) {
		return endianSwapIndex(x);
	} else {
		return x;
	}
}

/**
 * Read a unsigned from a file descriptor, optionally inverting
 * endianness.
 */
#ifdef BOWTIE_MM
template <typename index_t>
static inline index_t readIndex(int in, bool swap) {
	index_t x;
	if(read(in, (void *)&x, sizeof(index_t)) != sizeof(index_t)) {
		assert(false);
	}
	if(swap) {
		if(sizeof(index_t) == 8) {
			assert(false);
			return 0;
		} else if(sizeof(index_t) == 4) {
			return endianSwapU32(x);
		} else {
			assert_eq(sizeof(index_t), 2);
			return endianSwapU16(x);
		}
	} else {
		return x;
	}
}
#endif

/**
 * Read a unsigned from a FILE*, optionally inverting
 * endianness.
 */
template <typename index_t>
static inline index_t readIndex(FILE* in, bool swap) {
	index_t x;
	if(fread((void *)&x, 1, sizeof(index_t), in) != sizeof(index_t)) {
		assert(false);
	}
	if(swap) {
		if(sizeof(index_t) == 8) {
			assert(false);
			return 0;
		} else if(sizeof(index_t) == 4) {
			return endianSwapU32(x);
		} else {
			assert_eq(sizeof(index_t), 2);
			return endianSwapU16(x);
		}
	} else {
		return x;
	}
}


#endif /*WORD_IO_H_*/

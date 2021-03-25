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

/*
 * sse_wrap.h
 *
 * Routines to wrap Streaming SIMD Extensions (SSE) emmintrin.h
 * for an Intel x86 CPU and SIMD Everywhere (simde) for other CPUs.
 */

#ifndef SSE_WRAP_H_
#define SSE_WRAP_H_

#if defined(__aarch64__) || defined(__s390x__) || defined(__powerpc__)
#define SIMDE_ENABLE_NATIVE_ALIASES
#include "simde/x86/sse2.h"
#else
#include <emmintrin.h>
#endif

#endif /* SSE_WRAP_H_ */

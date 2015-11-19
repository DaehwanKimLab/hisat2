/*
 * Copyright 2015, Daehwan Kim <infphilo@gmail.com>
 *
 * This file is part of HISAT 2.
 *
 * Beast is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HISAT 2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with HISAT 2.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef HIERGBWT_COMMON_H_
#define HIERGBWT_COMMON_H_

// maximum size of a sequence represented by a local index
static const uint32_t local_index_size = (1 << 16) - (1 << 13);  // 1 << 5 is necessary for eftab index
static const uint32_t local_max_gbwt = (1 << 16) - (1 << 11);

// size of the overlapped sequence between the sequences represented by two consecutive local indexes
static const uint32_t local_index_overlap  = 1024;

// interval between two consecutive local indexes 
static const uint32_t local_index_interval = local_index_size - local_index_overlap;

// line rate in local indexes
static const int32_t local_lineRate_fm = 6;
static const int32_t local_lineRate_gfm = 7;

// how many rows are marked in a local index, every 2^<int>th row is marked
static const int32_t local_offRate = 3;

// the look table in a local index 4^<int> entries
static const int32_t local_ftabChars = 6;

#endif /*HIERGBWT_COMMON_H_*/

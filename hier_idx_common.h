/*
 * Copyright 2013, Daehwan Kim <infphilo@gmail.com>
 *
 * This file is part of Beast.  Beast is based on Bowtie 2.
 *
 * Beast is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Beast is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Beast.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef HIEREBWT_COMMON_H_
#define HIEREBWT_COMMON_H_

// maximum size of a sequence represented by a local index
static const uint32_t local_index_size     = (1 << 16) - (1 << 8);  // 1 << 5 is necessary for eftab index

// size of the overlapped sequence between the sequences represented by two consecutive local indexes
static const uint32_t local_index_overlap  = 1024;

// interval between two consecutive local indexes 
static const uint32_t local_index_interval = local_index_size - local_index_overlap;

// line rate in local indexes
static const int32_t local_lineRate = 6;

// how many rows are marked in a local index, every 2^<int>th row is marked
static const int32_t  local_offRate        = 3;

// the look table in a local index 4^<int> entries
static const int32_t  local_ftabChars      = 6;

#endif /*HIEREBWT_COMMON_H_*/

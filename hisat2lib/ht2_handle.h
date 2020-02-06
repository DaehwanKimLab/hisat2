/*
 * Copyright 2018, Chanhee Park <parkchanhee@gmail.com> and Daehwan Kim <infphilo@gmail.com>
 *
 * This file is part of HISAT 2.
 *
 * HISAT 2 is free software: you can redistribute it and/or modify
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

#ifndef __HT2_HANDLE_H__
#define __HT2_HANDLE_H__

#if 1
#define EXPORT __attribute__((visibility("default")))
#else
#define EXPORT 
#endif

typedef TIndexOffU index_t;
typedef uint16_t local_index_t;

struct ht2_handle {
    ALTDB<TIndexOffU>*      altdb;
    ALTDB<TIndexOffU>*      raltdb;
    RepeatDB<TIndexOffU>*   repeatdb;

    HGFM<TIndexOffU, local_index_t>* gfm;
    RFM<TIndexOffU> *rgfm;

    string tmp_str;

    string ht2_idx_name;

    struct ht2_options options; 
};

#endif /* __HT2_HANDLE_H__ */

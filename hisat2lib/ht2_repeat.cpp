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

#include <iostream>

#include "ds.h"
#include "repeat.h"
#include "rfm.h"

#include "ht2.h"
#include "ht2_handle.h"

EXPORT
void ht2_repeat_dump_repeatmap(ht2_handle_t handle)
{
    struct ht2_handle *hp = (struct ht2_handle *)handle;

    size_t localRFMSize = hp->rgfm->_localRFMs.size();

    cerr << "ht2lib: " << "LocalRFM size: " << localRFMSize << endl;
    cerr << "ht2lib: " << "Dump repeatMap" << endl;

    // repID -> 0
    const EList<pair<index_t, index_t> >& repeatMap = hp->repeatdb->repeatMap()[0];

    cerr << repeatMap.size() << endl;

    for(size_t i = 0; i < repeatMap.size(); i++) {
        cerr << repeatMap[i].first << ", " << repeatMap[i].second << endl;
    }

}


EXPORT
ht2_error_t ht2_repeat_expand(ht2_handle_t handle,
        const char *repeat_name,
        uint64_t repeat_pos,
        uint64_t repeat_len,
        struct ht2_repeat_expand_result **result_ptr)
{
    struct ht2_handle *hp = (struct ht2_handle *)handle;

    index_t rep_id = hp->rgfm->getLocalRFM_idx(repeat_name);

    TIndexOffU left = repeat_pos;
    TIndexOffU right = left + repeat_len;

    bool ret = hp->repeatdb->repeatExist(rep_id, left, right);
    if(!ret) {
        return HT2_ERR_NOT_REPEAT;
    }


    /* get coord */
    EList<pair<RepeatCoord<index_t>, RepeatCoord<index_t> > > positions;

    EList<index_t> snp_id_list;
    snp_id_list.clear();

    hp->repeatdb->getCoords(
            rep_id,
            left, right,
            snp_id_list,
            *(hp->raltdb),
            positions
            );

    /* build result */
    size_t result_size = sizeof(struct ht2_repeat_expand_result) + positions.size() * sizeof(struct ht2_position);

    struct ht2_repeat_expand_result *result = (struct ht2_repeat_expand_result *)malloc(result_size);

    result->count = positions.size();

    for(size_t i = 0; i < positions.size(); i++) {
        const RepeatCoord<index_t>& coord = positions[i].first;

        result->positions[i].chr_id = coord.tid;
        result->positions[i].pos = coord.toff;
        result->positions[i].direction = coord.fw ? 0 : 1; 
    }

    *(result_ptr) = result;
    return HT2_OK;
}


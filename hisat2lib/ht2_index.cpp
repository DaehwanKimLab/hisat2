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
const char* ht2_index_getrefnamebyid(ht2_handle_t handle, uint32_t chr_id)
{
    struct ht2_handle *hp = (struct ht2_handle *)handle;

    size_t refname_size = hp->gfm->_refnames.size();

    if(chr_id < (refname_size - 1)) {
        return hp->gfm->_refnames[chr_id].c_str();
    }

    return NULL;
}



EXPORT
ht2_error_t ht2_index_getrefnames(ht2_handle_t handle, struct ht2_index_getrefnames_result **result_ptr)
{
    struct ht2_handle *hp = (struct ht2_handle *)handle;

    size_t refname_size = hp->gfm->_refnames.size();
    if(refname_size == 0) {
        return HT2_ERR;
    }

    size_t result_hdr_size = sizeof(struct ht2_index_getrefnames_result) + sizeof(char *) * refname_size;
    size_t result_buf_size = 0;

    for(size_t i = 0; i < refname_size - 1; i++) {
        result_buf_size += strlen(hp->gfm->_refnames[i].c_str()) + 1;
    }

    void* ptr = malloc(result_hdr_size + result_buf_size);
    char* buf_ptr = (char *)ptr + result_hdr_size;

    memset(ptr, 0, result_hdr_size + result_buf_size);
    struct ht2_index_getrefnames_result* result = (struct ht2_index_getrefnames_result *)ptr;

    result->count = refname_size - 1;
    result->names[0] = buf_ptr;
    for(size_t i = 0; i < refname_size - 1; i++) {
        size_t rlen = strlen(hp->gfm->_refnames[i].c_str());
        strcpy(result->names[i], hp->gfm->_refnames[i].c_str());
        result->names[i + 1] = result->names[i] + rlen + 1;
    }

    (*result_ptr) = result;

    return HT2_OK;
}

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

using namespace std;

static const struct ht2_options ht2_default_options = {
    .offRate = -1,

    .useMm = false,
    .useShmem = false,
    .mmSweep = false,
    .noRefNames = false,
    .noSplicedAlignment = false, 
    .gVerbose = false, 
    .startVerbose = false,
    .sanityCheck = 0,

    .useHaplotype = false,
};

static void free_handle(struct ht2_handle *hp)
{
    if(hp->altdb) {
        delete hp->altdb;
    }

    if(hp->raltdb) {
        delete hp->raltdb;
    }

    if(hp->repeatdb) {
        delete hp->repeatdb;
    }

    if(hp->gfm) {
        delete hp->gfm;
    }

    if(hp->rgfm) {
        delete hp->rgfm;
    }

    delete hp;
}


static void init_handle(struct ht2_handle *hp)
{

    struct ht2_options *opt = &hp->options;

    hp->altdb = new ALTDB<index_t>();
    hp->repeatdb = new RepeatDB<TIndexOffU>();
    hp->raltdb = new ALTDB<index_t>();

    hp->gfm = new HGFM<TIndexOffU, local_index_t>(
            hp->ht2_idx_name,
            hp->altdb,
            NULL,
            NULL,
            -1,
            true,
            opt->offRate,
            0,
            opt->useMm,
            opt->useShmem,
            opt->mmSweep,
            !opt->noRefNames,
            true,
            true,
            true,
            !opt->noSplicedAlignment,
            opt->gVerbose,
            opt->startVerbose,
            false,
            opt->sanityCheck,
            opt->useHaplotype);


    // Load the other half of the index into memory
    assert(!hp->gfm->isInMemory());

    hp->gfm->loadIntoMemory(
            -1,
            true,
            true,
            true,
            !opt->noRefNames,
            opt->startVerbose);

    hp->rgfm = new RFM<TIndexOffU>(
            hp->ht2_idx_name + ".rep",
            hp->raltdb,
            hp->repeatdb,
            NULL,
            -1,
            true,
            opt->offRate,
            0,
            opt->useMm,
            opt->useShmem,
            opt->mmSweep,
            !opt->noRefNames,
            true,
            true,
            true,
            !opt->noSplicedAlignment,
            opt->gVerbose,
            opt->startVerbose,
            false,
            opt->sanityCheck,
            false);


    assert(!hp->rgfm->isInMemory());
    hp->rgfm->loadIntoMemory(
            -1,
            true,
            true,
            true,
            !opt->noRefNames,
            opt->startVerbose);

    hp->repeatdb->construct(hp->gfm->rstarts(), hp->gfm->nFrag());
}


ht2_handle_t ht2_init(const char *name, ht2_option_t *options)
{
    cerr << "ht2_init" << endl;

    struct ht2_handle *handle = new ht2_handle;

    handle->ht2_idx_name = name;
    if(options) {
        memcpy(&handle->options, options, sizeof(struct ht2_options));
    } else {
        memcpy(&handle->options, &ht2_default_options, sizeof(struct ht2_options));
    }

    // Init
    init_handle(handle);

    handle->tmp_str = name;

    return (ht2_handle_t)handle;
}


void ht2_close(ht2_handle_t handle)
{
    struct ht2_handle *hp = (struct ht2_handle *)handle;
    if(hp == NULL) {
        return;
    }

    cerr << __func__ << ": " << hp->tmp_str << endl;

    //delete str;

    // free memory

    free_handle(hp);
}



ht2_error_t ht2_init_options(ht2_option_t *options)
{
    if(options == NULL) {
        return HT2_ERR;
    }

    memcpy(options, &ht2_default_options, sizeof(ht2_default_options));

    return HT2_OK;
}


void ht2_test_1(ht2_handle_t handle)
{
    struct ht2_handle *hp = (struct ht2_handle *)handle;

    size_t refname_size = hp->gfm->_refnames.size();
    cerr << "ht2lib: " << "gfm refnames: " << refname_size << endl;
    for(size_t i = 0; i < refname_size; i++) {
        cerr << "ht2lib: " << " " << i << " -> " << hp->gfm->_refnames[i] << endl;
    }
}

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

#ifndef __HT2_H__
#define __HT2_H__


#ifdef __cplusplus
extern "C"
{
#endif


typedef int ht2_error_t;

enum {
    HT2_OK              = 0,

    HT2_ERR             = -1,
    HT2_ERR_NOT_REPEAT  = -2,
};
#define HT2_RET_OK(x)   ((x) == HT2_OK)

typedef void* ht2_handle_t;

struct ht2_options {
    int offRate;

    bool useMm;
    bool useShmem;
    bool mmSweep;
    bool noRefNames;
    bool noSplicedAlignment;
    int gVerbose;
    bool startVerbose;
    int sanityCheck;
    
    bool useHaplotype;
};

typedef struct ht2_options ht2_option_t;



/**************************************************************************
 *
 * Initialize APIs
 *
 **************************************************************************/

ht2_handle_t ht2_init(const char *name, ht2_option_t *options);
void ht2_close(ht2_handle_t);

ht2_error_t ht2_init_options(ht2_option_t *options);


/**************************************************************************
 *
 * Index APIs
 *
 **************************************************************************/

const char* ht2_index_getrefnamebyid(ht2_handle_t handle, uint32_t chr_id);  

struct ht2_index_getrefnames_result {
    int count;
    char* names[0];
};

/**
 * @brief 
 *
 * @param handle
 * @param result_ptr        pointer to result. Caller must relese memory by free().
 *
 * @return 
 */
ht2_error_t ht2_index_getrefnames(ht2_handle_t handle, struct ht2_index_getrefnames_result **result_ptr);


/**************************************************************************
 *
 * Repeat APIs
 *
 **************************************************************************/

struct ht2_position {
    uint32_t chr_id;
    int direction; /* 0 - forward, 1 - reverse */
    uint64_t pos;
};

struct ht2_repeat_expand_result {
    int count;
    struct ht2_position positions[0];
};


/**
 * @brief 
 *
 * @param handle
 * @param repeat_name
 * @param repeat_pos        repeat position on repeat sequence(1-based)
 * @param repeat_len
 * @param result_ptr        pointer to result. caller must release memory by free(). 
 *                          ex) free(result_ptr);
 *
 * @return 
 */
ht2_error_t ht2_repeat_expand(ht2_handle_t handle, 
        const char *repeat_name, 
        uint64_t repeat_pos, 
        uint64_t repeat_len, 
        struct ht2_repeat_expand_result **result_ptr);

/**************************************************************************
 *
 * Alignment APIs
 *
 **************************************************************************/
/* TODO */


/**************************************************************************
 *
 * ETC APIs
 *
 **************************************************************************/
/* TODO */
void ht2_test_1(ht2_handle_t);
void ht2_repeat_dump_repeatmap(ht2_handle_t handle);


#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* __HT2_H__ */

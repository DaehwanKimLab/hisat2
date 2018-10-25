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
#include <stdarg.h>
#include <Python.h>

#include "ht2.h"

#define HT2_HANDLE_ID "handle"

#ifdef DEBUG
#define DEBUGLOG(fmt, ...) do { fprintf(stderr, "%s:%d:%s(): " fmt, __FILE__, __LINE__, __func__, ##__VA_ARGS__);  } while(0) 
#else
#define DEBUGLOG(fmt, ...) 
#endif

static ht2_handle_t get_handle(PyObject *cap)
{
	return PyCapsule_GetPointer(cap, HT2_HANDLE_ID);
}

static PyObject *conv_refnames_result(struct ht2_index_getrefnames_result *result)
{
    PyObject *refnames = NULL;
    size_t i = 0;

    if(result == NULL) {
        return NULL;
    }

    refnames = PyList_New(result->count);
    for(i = 0; i < result->count; i++) {
        PyObject *str = PyString_FromString(result->names[i]);
        PyList_SetItem(refnames, i, str);
    }

    return refnames;
}

static PyObject *conv_repeat_expand_result(struct ht2_repeat_expand_result *result)
{
    PyObject *positions = NULL;
    size_t i = 0;

    if(result == NULL) {
        return NULL;
    }

    positions = PyList_New(result->count);
    for(i = 0; i < result->count; i++) {
        struct ht2_position *htpos = &result->positions[i];

        PyList_SetItem(positions, i, 
                Py_BuildValue("(III)", htpos->chr_id, htpos->direction, htpos->pos)
                );
    }

    return positions;
}


static PyObject *conv_ht2opt(ht2_option_t *opts)
{
	PyObject *py_opt = NULL;

	py_opt = PyDict_New();
	if(py_opt == NULL) {
		return NULL;
	}
#define HT2_OPT_BUILD(_pobj, _popt, _name, _type) \
	do {\
		if(PyDict_SetItemString((_pobj), #_name, Py_BuildValue((_type), (_popt)->_name)) < 0) {\
			DEBUGLOG("Can't set item %s\n", #_name);\
		} \
	} while(0)


	HT2_OPT_BUILD(py_opt, opts, offRate, "i");
	HT2_OPT_BUILD(py_opt, opts, useMm, "i");
	HT2_OPT_BUILD(py_opt, opts, useShmem, "i");
	HT2_OPT_BUILD(py_opt, opts, mmSweep, "i");
	HT2_OPT_BUILD(py_opt, opts, noRefNames, "i");
	HT2_OPT_BUILD(py_opt, opts, noSplicedAlignment, "i");
	HT2_OPT_BUILD(py_opt, opts, gVerbose, "i");
	HT2_OPT_BUILD(py_opt, opts, startVerbose, "i");
	HT2_OPT_BUILD(py_opt, opts, sanityCheck, "i");
	HT2_OPT_BUILD(py_opt, opts, useHaplotype, "i");

	return py_opt;
}

static void update_ht2_options(ht2_option_t *ht2opt, PyObject *py_opt)
{
#define HT2_OPT_UPDATE(_pobj, _ht2opt, _name) \
	do {\
		PyObject *p;\
		if((p = PyDict_GetItemString((_pobj), #_name)) != NULL) { \
			(_ht2opt)->_name = PyInt_AsLong(p); \
			DEBUGLOG(#_name " %d\n", (ht2opt)->_name); \
			if(PyErr_Occurred() != NULL) { \
				DEBUGLOG("Error Occurred"); \
			}\
		}\
	} while (0)

	HT2_OPT_UPDATE(py_opt, ht2opt, offRate);
	HT2_OPT_UPDATE(py_opt, ht2opt, useMm);
	HT2_OPT_UPDATE(py_opt, ht2opt, mmSweep);
	HT2_OPT_UPDATE(py_opt, ht2opt, noRefNames);
	HT2_OPT_UPDATE(py_opt, ht2opt, noSplicedAlignment);
	HT2_OPT_UPDATE(py_opt, ht2opt, gVerbose);
	HT2_OPT_UPDATE(py_opt, ht2opt, startVerbose);
	HT2_OPT_UPDATE(py_opt, ht2opt, sanityCheck);
	HT2_OPT_UPDATE(py_opt, ht2opt, useHaplotype);

}

static PyObject *ht2py_get_options(PyObject *self, PyObject *args)
{
	ht2_option_t ht2opt;

	ht2_init_options(&ht2opt);


	/* convert ht2_option_t to PyObject(map) */

	PyObject *pobj = conv_ht2opt(&ht2opt);

	return pobj;
}

static PyObject *ht2py_init(PyObject *self, PyObject *args)
{
	ht2_handle_t handle;
	PyObject *popt = NULL;
	char *name = NULL;

	if(!PyArg_ParseTuple(args, "sO", &name, &popt)) {
		return NULL;
	}

	DEBUGLOG("name %s\n", name);
	DEBUGLOG("popt %p\n", popt);

	if(!PyDict_CheckExact(popt)) {
		// TODO
		// exception
		DEBUGLOG("Invalid data type\n");
		return NULL;
	}

	ht2_option_t ht2opt;
	ht2_init_options(&ht2opt);
	update_ht2_options(&ht2opt, popt);

	handle = ht2_init(name, &ht2opt);

	DEBUGLOG("handle %p\n", handle);

	PyObject *cap = PyCapsule_New(handle, HT2_HANDLE_ID, NULL);

	return cap;
}

static PyObject *ht2py_close(PyObject *self, PyObject *args)
{
	ht2_handle_t handle;
    PyObject *cap;

    // Parse Args
    // ht2py.close(handle)
    //
	if(!PyArg_ParseTuple(args, "O", &cap)) {
		DEBUGLOG("Can't parse args\n");
		return NULL;
	}

	handle = get_handle(cap);
	if(handle == NULL) {
		DEBUGLOG("Can't get handle\n");
		return NULL;
	}

	DEBUGLOG("handle %p\n", handle);

	ht2_close(handle);

	Py_RETURN_NONE;
}


static PyObject *ht2py_index_getrefnamebyid(PyObject *self, PyObject *args)
{
    PyObject *cap;
    uint32_t chr_id;

    // ht2py.index_getrefnamebyid(handle, chr_id)

    if(!PyArg_ParseTuple(args, "Oi", &cap, &chr_id)) {
		DEBUGLOG("Can't parse args\n");
		return NULL;
    }

	ht2_handle_t handle = get_handle(cap);
	if(handle == NULL) {
		DEBUGLOG("Can't get handle\n");
		return NULL;
	}

    const char *refname = ht2_index_getrefnamebyid(handle, chr_id);

    if(refname == NULL) {
		DEBUGLOG("Can't get refname(%u)\n", chr_id);
        return Py_BuildValue("s", "");
    }

    return Py_BuildValue("s", refname);
}

static PyObject *ht2py_index_getrefnames(PyObject *self, PyObject *args)
{
    ht2_handle_t handle;
    PyObject *cap;

    // Parse Args
    // ht2py.index_getrefnames(handle)
    if(!PyArg_ParseTuple(args, "O", &cap)) {
		DEBUGLOG("Can't parse args\n");
        return NULL;
    }

    handle = get_handle(cap);
    if(handle == NULL) {
		DEBUGLOG("Can't get handle\n");
        return NULL;
    }


    struct ht2_index_getrefnames_result *result = NULL;
    ht2_error_t ret = ht2_index_getrefnames(handle, &result);

    PyObject *refnames = NULL;
    
    if(ret == HT2_OK) {
        /* Build List of names */
        refnames = conv_refnames_result(result);
        free(result);
    } else {
        refnames = PyList_New(0);
    }

    return refnames;
}

static PyObject *ht2py_repeat_expand(PyObject *self, PyObject *args)
{
    PyObject *cap;
    char *name = NULL;
    uint64_t rpos = 0;
    uint64_t rlen = 0;

    // Parse Args
    // ht2py.repeat_expand(handle, 'repeat_name', repeat_pos, repeat_len)
    if(!PyArg_ParseTuple(args, "OsLL", &cap, &name, &rpos, &rlen)) {
		DEBUGLOG("Can't parse args\n");
        return NULL;
    }

    //fprintf(stderr, "%s, %lu, %lu\n", name, rpos, rlen);

    ht2_handle_t handle = get_handle(cap);
    if(handle == NULL) {
		DEBUGLOG("Can't get handle\n");
        return NULL;
    }


    struct ht2_repeat_expand_result *result = NULL;
    ht2_error_t ret = ht2_repeat_expand(handle, name, rpos, rlen, &result);

    PyObject *positions = NULL;
    if(ret == HT2_OK) {
        /* Build list of position */
        positions = conv_repeat_expand_result(result);
        free(result);
    } else {
		DEBUGLOG("error %d, %s, %lu, %lu\n", ret,
				name, rpos, rlen);
        positions = PyList_New(0);
    }

    return positions;
}


static PyMethodDef myMethods[] = {
	/* Initialize APIs */
	{"get_options", ht2py_get_options, METH_NOARGS, "Get default options"},
	{"init", ht2py_init, METH_VARARGS, "Initialize HT2Lib handle"},
	{"close", ht2py_close, METH_VARARGS, "Release HT2Lib handle"},

	/* Index APIs */
	{"index_getrefnamebyid", ht2py_index_getrefnamebyid, METH_VARARGS, "Get reference name"},
	{"index_getrefnames", ht2py_index_getrefnames, METH_VARARGS, "Get all reference names"},

	/* Repeat APIs */
	{"repeat_expand", ht2py_repeat_expand, METH_VARARGS, "Find reference positions"},

	/* */
	{NULL, NULL, 0, NULL}
};


PyMODINIT_FUNC
initht2py(void)
{
	(void)Py_InitModule("ht2py", myMethods);
}

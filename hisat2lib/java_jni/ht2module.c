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

#include <jni.h>
#include <stdio.h>
#include "ht2.h"

#include "HT2Module.h"

#define CLASSPATH_INTEGER		"java/lang/Integer"
#define CLASSPATH_HASHMAP		"java/util/HashMap" 
#define CLASSPATH_ARRAYLIST		"java/util/ArrayList" 
#define CLASSPATH_HT2POSITION	"HT2Module$HT2Position"

#ifdef DEBUG
#define DEBUGLOG(fmt, ...) do { fprintf(stderr, "%s:%d:%s(): " fmt, __FILE__, __LINE__, __func__, ##__VA_ARGS__);  } while(0) 
#else
#define DEBUGLOG(fmt, ...) 
#endif

static jint JNI_VERSION = JNI_VERSION_1_2;

static jclass classInteger;
static jclass classHashMap;
static jclass classArrayList;
static jclass classHT2Position;

static jmethodID mthIntegerInit;
static jmethodID mthIntegerIntValue;

static jmethodID mthHashMapInit;
static jmethodID mthHashMapPut;
static jmethodID mthHashMapGet;

static jmethodID mthArrayListInit;
static jmethodID mthArrayListEnsureCapacity;
static jmethodID mthArrayListAdd;

static jmethodID mthHT2PositionInit;

static jobject NewHashMap(JNIEnv *env)
{
	jobject obj = (*env)->NewObject(env, classHashMap, mthHashMapInit);
	return obj;
}

static jobject NewInteger(JNIEnv *env, int value)
{
	jobject obj = (*env)->NewObject(env, classInteger, mthIntegerInit, value);
	return obj;
}

static jobject NewArrayList(JNIEnv *env)
{
	jobject obj = (*env)->NewObject(env, classArrayList, mthArrayListInit);
	return obj;
}

static jobject NewHT2Position(JNIEnv *env, struct ht2_position *htpos)
{
	jobject obj = (*env)->NewObject(env, classHT2Position, mthHT2PositionInit, 
			htpos->chr_id, htpos->direction, htpos->pos);
	return obj;
}

static int GetInteger(JNIEnv *env, jobject jobjInt, int *val)
{
	int value = (*env)->CallIntMethod(env, jobjInt, mthIntegerIntValue);
	*val = value;
	return 0;
}

static void hashmap_put(JNIEnv *env, jobject jobjMap,
		const char *key_str, const int val)
{
	jstring keyObj;
	jobject valueObj;

	keyObj = (*env)->NewStringUTF(env, key_str); 
	valueObj = NewInteger(env, val);

	(*env)->CallObjectMethod(env, jobjMap, mthHashMapPut, keyObj, valueObj);

	(*env)->DeleteLocalRef(env, keyObj);
	(*env)->DeleteLocalRef(env, valueObj);
}

static int hashmap_get(JNIEnv *env, jobject jobjMap, 
		const char *key_str, int *val)
{
	jstring keyObj;
	jobject valueObj;

	keyObj = (*env)->NewStringUTF(env, key_str);
	valueObj = (*env)->CallObjectMethod(env, jobjMap, mthHashMapGet, keyObj);

	if(valueObj == NULL) {
		DEBUGLOG("Can't find key: %s\n", key_str);
		(*env)->DeleteLocalRef(env, keyObj);
		return -1;
	}

	// Integer -> int
	int value = 0;	
	if(GetInteger(env, valueObj, &value) < 0) {
		DEBUGLOG("Can't get Integer\n");
		(*env)->DeleteLocalRef(env, keyObj);
		(*env)->DeleteLocalRef(env, valueObj);
		return -1;
	}

	(*env)->DeleteLocalRef(env, keyObj);
	(*env)->DeleteLocalRef(env, valueObj);

	*val = value;

	return 0;
}

static jobject conv_ht2option(JNIEnv *env, ht2_option_t *opts)
{
	jobject hashMap = NewHashMap(env);

#define HT2_OPT_BUILD(_name) \
	hashmap_put(env, hashMap, #_name , opts->_name)

	HT2_OPT_BUILD(offRate);
	HT2_OPT_BUILD(useMm);
	HT2_OPT_BUILD(useShmem);
	HT2_OPT_BUILD(mmSweep);
	HT2_OPT_BUILD(noRefNames);
	HT2_OPT_BUILD(noSplicedAlignment);
	HT2_OPT_BUILD(gVerbose);
	HT2_OPT_BUILD(startVerbose);
	HT2_OPT_BUILD(sanityCheck);
	HT2_OPT_BUILD(useHaplotype);

	return hashMap;
}

static void update_ht2option(JNIEnv *env, ht2_option_t *opts, jobject jmapObject)
{
#define HT2_OPT_UPDATE(_name) \
	do { \
		int value = 0; \
		if(hashmap_get(env, jmapObject, #_name, &value) == 0) { \
			DEBUGLOG("Using %s, %d\n", #_name, value); \
			opts->_name = value; \
		} else { \
			DEBUGLOG("Using default %s, %d\n", #_name, opts->_name); \
		}\
	} while(0)

	HT2_OPT_UPDATE(offRate);
	HT2_OPT_UPDATE(useMm);
	HT2_OPT_UPDATE(useShmem);
	HT2_OPT_UPDATE(mmSweep);
	HT2_OPT_UPDATE(noRefNames);
	HT2_OPT_UPDATE(noSplicedAlignment);
	HT2_OPT_UPDATE(gVerbose);
	HT2_OPT_UPDATE(startVerbose);
	HT2_OPT_UPDATE(sanityCheck);
	HT2_OPT_UPDATE(useHaplotype);
}


static void conv_refnames_result(JNIEnv *env, 
		struct ht2_index_getrefnames_result *result, jobject jobjList)
{
	size_t i;

	// Resize
	DEBUGLOG("count: %d\n", result->count);
	(*env)->CallVoidMethod(env, jobjList, mthArrayListEnsureCapacity, result->count);

	for(i = 0; i < result->count; i++) {
		DEBUGLOG("names[%d] %s\n", i, result->names[i]);

		jobject elem = (*env)->NewStringUTF(env, result->names[i]);
		(*env)->CallObjectMethod(env, jobjList, mthArrayListAdd, elem);
		(*env)->DeleteLocalRef(env, elem);
	}
}

static void conv_repeat_expand_result(JNIEnv *env,
		struct ht2_repeat_expand_result *result,
		jobject jobjList)
{
	size_t i;

	// Resize
	DEBUGLOG("count: %d\n", result->count);
	(*env)->CallVoidMethod(env, jobjList, mthArrayListEnsureCapacity, result->count);

	// Add Items
	for(i = 0; i < result->count; i++) {
        struct ht2_position *htpos = &result->positions[i];

		DEBUGLOG("position[%d]: %u, %d, %lu\n", i, htpos->chr_id, htpos->direction, htpos->pos);

		jobject elem = NewHT2Position(env, htpos);
		(*env)->CallObjectMethod(env, jobjList, mthArrayListAdd, elem);
		(*env)->DeleteLocalRef(env, elem);
	}
}

JNIEXPORT jlong JNICALL Java_HT2Module_init(JNIEnv *env,
		jobject thisObj, jstring indexNameObj, jobject optionMap)
{
	DEBUGLOG("Init\n");

	const char *index_path = (*env)->GetStringUTFChars(env, indexNameObj, NULL);
	if(index_path == NULL) {
		DEBUGLOG("Can't get string\n");
		return 0;
	}

	DEBUGLOG("Index Path: %s\n", index_path);

	ht2_option_t ht2opt;
	ht2_init_options(&ht2opt);
	if(optionMap != NULL) {
		// update options
		update_ht2option(env, &ht2opt, optionMap);
	}

	ht2_handle_t handle = ht2_init(index_path, &ht2opt);

	DEBUGLOG("Initailzed %p\n", handle);

	(*env)->ReleaseStringUTFChars(env, indexNameObj, index_path);

	return (jlong)handle;
}


JNIEXPORT void JNICALL Java_HT2Module_close(JNIEnv *env,
		jobject thisObj, jlong handlePtr)
{
	DEBUGLOG("Closing\n");

	ht2_handle_t handle = (ht2_handle_t)handlePtr;

	DEBUGLOG("Received handle: %p\n", handle);

	ht2_close(handle);
}

JNIEXPORT jobject JNICALL Java_HT2Module_get_1options(JNIEnv *env,
		jobject thisObj)
{

	DEBUGLOG("get_option\n");

	// Get Default options 
	ht2_option_t ht2opt;
	ht2_init_options(&ht2opt);

	// convert ht2opt to java hashmap
	jobject hashobj = conv_ht2option(env, &ht2opt);

	return hashobj;
}

JNIEXPORT jstring JNICALL Java_HT2Module_index_1getrefnamebyid(JNIEnv *env,
		jobject thisObj, jlong handlePtr, jint chr_id)
{
	DEBUGLOG("index_getrefnamebyid\n");
	ht2_handle_t handle = (ht2_handle_t)handlePtr;
	if(handle == NULL) {
		DEBUGLOG("Invalid handle\n");
		return NULL;
	}

	const char *refname = ht2_index_getrefnamebyid(handle, chr_id);
	if(refname == NULL) {
		DEBUGLOG("Can't get refname(%d)\n", chr_id);
		return NULL;
	}

	return (*env)->NewStringUTF(env, refname);
}

JNIEXPORT jobject JNICALL Java_HT2Module_index_1getrefnames(JNIEnv *env,
		jobject thisObj, jlong handlePtr)
{
	DEBUGLOG("index_getrefnames\n");

	ht2_handle_t handle = (ht2_handle_t)handlePtr;
	if(handle == NULL) {
		DEBUGLOG("Invalid handle\n");
		return NULL;
	}

	struct ht2_index_getrefnames_result *result = NULL;
	ht2_error_t ret = ht2_index_getrefnames(handle, &result);

	jobject refnames = NewArrayList(env);

	if(ret == HT2_OK) {
		DEBUGLOG("refnames size: %d\n", result->count);
		conv_refnames_result(env, result, refnames);
		free(result);
	} else {
		DEBUGLOG("Can't get refnames\n");
	}

	return refnames;
}

JNIEXPORT jobject JNICALL Java_HT2Module_repeat_1expand(JNIEnv *env,
		jobject thisObj, jlong handlePtr, jstring nameObj, jlong rpos, jlong rlen)
{
	DEBUGLOG("repeat_expand\n");

	ht2_handle_t handle = (ht2_handle_t)handlePtr;
	if(handle == NULL) {
		DEBUGLOG("Invalid handle\n");
		return NULL;
	}

	const char *name = (*env)->GetStringUTFChars(env, nameObj, NULL);
	if(name == NULL) {
		DEBUGLOG("Can't get string\n");
		return NULL;
	}

	DEBUGLOG("Repeat Expand: %s, %lu, %lu\n", name, rpos, rlen);


	struct ht2_repeat_expand_result *result = NULL;
	ht2_error_t ret = ht2_repeat_expand(handle, name, rpos, rlen, &result);

	jobject positions = NewArrayList(env);

	if(ret == HT2_OK) {
		DEBUGLOG("expand position size: %d\n", result->count);
		conv_repeat_expand_result(env, result, positions);
		free(result);
	} else {
		DEBUGLOG("Can't expand repeat\n");
	}

	(*env)->ReleaseStringUTFChars(env, nameObj, name);

	return positions;
}

jint JNI_OnLoad(JavaVM *vm, void *reserved)
{
	DEBUGLOG("JNI Loaded\n");

	JNIEnv *env;

	if((*vm)->GetEnv(vm, (void **)&env, JNI_VERSION) != JNI_OK) {
		return JNI_ERR;
	}

	// Load Class
#define LOAD_CLASS(_globalRef, _clazz) \
	do { \
		jclass tmpClassRef = (*env)->FindClass(env, _clazz); \
		if(tmpClassRef == NULL) { \
			DEBUGLOG("Can't find class %s\n", _clazz); \
			return JNI_ERR; \
		} \
		_globalRef = (jclass)(*env)->NewGlobalRef(env, tmpClassRef); \
		(*env)->DeleteLocalRef(env, tmpClassRef); \
	} while(0)

	LOAD_CLASS(classInteger, CLASSPATH_INTEGER);
	LOAD_CLASS(classHashMap, CLASSPATH_HASHMAP);
	LOAD_CLASS(classArrayList, CLASSPATH_ARRAYLIST);
	LOAD_CLASS(classHT2Position, CLASSPATH_HT2POSITION);


	// Load Method
	
	// Integer
	mthIntegerInit = (*env)->GetMethodID(env, classInteger, "<init>", "(I)V");
	mthIntegerIntValue = (*env)->GetMethodID(env, classInteger, "intValue", "()I");

	// HashMap
	mthHashMapInit = (*env)->GetMethodID(env, classHashMap, "<init>", "()V");
	mthHashMapPut = (*env)->GetMethodID(env, classHashMap, "put",
			"(Ljava/lang/Object;Ljava/lang/Object;)Ljava/lang/Object;");
	mthHashMapGet = (*env)->GetMethodID(env, classHashMap, "get",
			"(Ljava/lang/Object;)Ljava/lang/Object;");

	// ArrayList
	mthArrayListInit = (*env)->GetMethodID(env, classArrayList, "<init>", "()V");
	mthArrayListEnsureCapacity = (*env)->GetMethodID(env, classArrayList, "ensureCapacity", "(I)V");
	mthArrayListAdd = (*env)->GetMethodID(env, classArrayList, "add", 
			"(Ljava/lang/Object;)Z");

	// HT2Position
	mthHT2PositionInit = (*env)->GetMethodID(env, classHT2Position, "<init>", "(III)V");

	return JNI_VERSION;
}


void JNI_OnUnload(JavaVM *vm, void *reserved)
{
	DEBUGLOG("JNI Unload\n");
	JNIEnv *env;
	
	(*vm)->GetEnv(vm, (void **)&env, JNI_VERSION);

	(*env)->DeleteGlobalRef(env, classInteger);
	(*env)->DeleteGlobalRef(env, classHashMap);
	(*env)->DeleteGlobalRef(env, classArrayList);
	(*env)->DeleteGlobalRef(env, classHT2Position);

}


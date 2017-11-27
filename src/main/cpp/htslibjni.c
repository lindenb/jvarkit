#include <jni.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <htslib/kstring.h>
#include <htslib/hts.h>
#include "htslibjni.h"

static void* safeMalloc(size_t size) {
	void* ptr = malloc(size);
	if(ptr==NULL) {
		fprintf(stderr,"cannot alloc memory.");
		exit(EXIT_FAILURE);
		}
	memset(ptr,0,size);
	return ptr;
	}

#define C_NAME(STR) _c_ ## STR
#define STR_TO_C(javaString) const char* C_NAME(javaString) = (*env)->GetStringUTFChars(env, javaString, 0)

#define STR_RELEASE(nativeString) (*env)->ReleaseStringUTFChars(env, javaString, C_NAME(javaString))
	
#define CLASS(NAME)  Java_com_github_lindenb_jvarkit_htslib_HtsFile_ ##  NAME 
jlong CLASS(_1open) (JNIEnv* env, jclass c, jstring filename, jstring mode) {
	jlong f=0L;
	STR_TO_C(filename);
	/*STR_TO_C(mode);
	htsFile* f= hts_open(C_VAR(filename),C_VAR(mode));
	STR_RELEASE(filename);
	STR_RELEASE(mode);*/
	return (jlong)f;
	}


void CLASS(_1close) (JNIEnv* env, jclass c, jlong ptr) {
	hts_close((htsFile*)ptr);
	}

#undef CLASS

#define CLASS(NAME)  Java_com_github_lindenb_jvarkit_htslib_KString_ ##  NAME 

jlong CLASS(_1create) (JNIEnv *env, jclass clazz) {
	return (long)safeMalloc(sizeof(kstring_t));
	}

jbyte CLASS(_1at)(JNIEnv* env, jclass clazz, jlong ptr, jint idx) {
	return ks_str((kstring_t*)ptr)[idx];
	}

jint CLASS(_len)(JNIEnv* env, jclass clazz, jlong ptr) {
	return (jint)ks_len((kstring_t*)ptr);
	}

void CLASS(_release)(JNIEnv* env, jclass clazz, jlong ptr) {
	kstring_t* ks = (kstring_t*)ptr;
	free(ks->s);
	free(ks);
	}


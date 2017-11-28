#include <jni.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <htslib/kstring.h>
#include <htslib/vcf.h>
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

#define C_VAR(S) _c_ ## S

#define STR_J2C(S) const char* C_VAR(S) = ( const char *) (*env)->GetStringUTFChars(env, S, NULL)
#define FREE_STR(S) (*env)->ReleaseStringUTFChars(env,S,C_VAR(S))

#define CLASS(NAME)  Java_com_github_lindenb_jvarkit_htslib_HtsFile_ ##  NAME 


jlong CLASS(_1open) (JNIEnv *env, jclass clazz, jstring filename, jstring mode) {
	STR_J2C(filename);
	STR_J2C(mode);
	htsFile* f= hts_open(C_VAR(filename),C_VAR(mode));
	FREE_STR(filename);
	FREE_STR(mode);
	return (jlong)f;
	}


void CLASS(_1close) (JNIEnv* env, jclass c, jlong ptr) {
	if(ptr<=0L) return;
	hts_close((htsFile*)ptr);
	}

#undef CLASS

#define CLASS(NAME)  Java_com_github_lindenb_jvarkit_htslib_KString_ ##  NAME 

jlong CLASS(_1create) (JNIEnv* env, jclass clazz) {
	return (jlong)safeMalloc(sizeof(kstring_t));
	}

jbyte CLASS(_1at)(JNIEnv* env, jclass clazz, jlong ptr, jint idx) {
	return (jbyte)ks_str((kstring_t*)ptr)[idx];
	}

jint CLASS(_len)(JNIEnv* env, jclass clazz, jlong ptr) {
	return (jint)ks_len((kstring_t*)ptr);
	}

void CLASS(_release)(JNIEnv* env, jclass clazz, jlong ptr) {
	kstring_t* ks = (kstring_t*)ptr;
	if(ptr==0L) return;
	free(ks->s);
	free(ks);
	}
#undef CLASS


#define CLASS(NAME)  Java_com_github_lindenb_jvarkit_htslib_Bcf1_ ##  NAME 

jlong CLASS(_1create) (JNIEnv* env, jclass clazz) {
	return (jlong)bcf_init();
	}

jboolean CLASS(_1is_1snp)(JNIEnv* env, jclass clazz, jlong adrs) {
	return (jboolean)bcf_is_snp((bcf1_t*)adrs);
	}

void CLASS(_1destroy)(JNIEnv* env, jclass clazz, jlong adrs) {
	bcf1_t* p = (bcf1_t*)adrs;
	if(p==0L) return;
	bcf_destroy(p);
	}
#undef CLASS



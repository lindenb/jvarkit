

#include <stdio.h>
#include <unistd.h>
#include <zlib.h>
#include <getopt.h>
#include "htslib/kseq.h"
#include <htslib/hts.h>

#include "com_github_lindenb_jvarkit_htslib_HtsLib.h"

KSEQ_DECLARE(gzFile)

#define QUALIFIEDMETHOD(fun) Java_com_github_lindenb_jvarkit_htslib_HtsLib_##fun

jlong QUALIFIEDMETHOD(kseq_1init_1file)(JNIEnv *env, jclass clazz, jstring filename) {
    const char* fn = (*env)->GetStringUTFChars(env,(jstring) filename, NULL);
    gzFile fp = strcmp(fn, "-") ? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
    if (fp == 0) {
        fprintf(stderr, "Cannot open %s.\n", fn);
        (*env)->ReleaseStringUTFChars(env,filename, fn);
        return (jlong)0L;
        }
    kseq_t* obj = kseq_init(fp);
    (*env)->ReleaseStringUTFChars(env,filename, fn);
    return (jlong)obj;
    }


void QUALIFIEDMETHOD(kseq_1destroy)(JNIEnv *env, jclass clazz,jlong ptr) {
  if(ptr!=0L)  kseq_destroy(( kseq_t*)ptr);
}

jint QUALIFIEDMETHOD(kseq_1read4y)(JNIEnv *env, jclass clazz,jlong ptr,jobjectArray array) {
  if(ptr==0L) return -1;
  kseq_t* kseq = ( kseq_t*)ptr;
  int rez = kseq_read(kseq);
  if(rez<0) return rez;
  jobject s =  (*env)->NewStringUTF(env,kseq->name.s);
  (*env)->SetObjectArrayElement(env,array, 0,s);
  s =  (*env)->NewStringUTF(env,kseq->seq.s);
  (*env)->SetObjectArrayElement(env,array, 1,s);
  s =  (*env)->NewStringUTF(env,kseq->comment.s);
  (*env)->SetObjectArrayElement(env,array, 2,s);
  s =  (*env)->NewStringUTF(env,kseq->qual.s);
  (*env)->SetObjectArrayElement(env,array, 3,s);
  return rez;
}


jint QUALIFIEDMETHOD(kseq_1read4y)(JNIEnv *env, jclass clazz,jlong ptr,jobjectArray array) {
  const char* fn = (*env)->GetStringUTFChars(env,(jstring) filename, NULL);
  htsFile *in = hts_open(fn,"r");
  if(in==NULL) {
    ERROR("Cannot open input vcf %s. (%s)",fn,strerror(errno));
    return -1;    
    }
 bcf_hdr_t *header = bcf_hdr_read(in);
 if(header==NULL) {
    ERROR("Cannot read header for input vcf %s.",fn);
    return -1;    
    }


 bcf1_t* bcf = bcf_init();
 if(bcf==NULL) {
    ERROR("Out of memory.");
    return -1;    
    }
  (*env)->ReleaseStringUTFChars(env,filename, fn);
  (*env)->SetLongArrayElement(env,array,0,(jlong)in);
  (*env)->SetLongArrayElement(env,array,1,(jlong)header);
  (*env)->SetLongArrayElement(env,array,2,(jlong)bcf);
  return 0;
  }

jint QUALIFIEDMETHOD(kseq_1read4y)(JNIEnv *env, jclass clazz,jlong ptr,jobjectArray array) {

htsFile *in = (htsFile*) (*env)->GetLongArrayElement(env,array,0);
if(in!=NULL) hts_close(in);

bcf_hdr_t *header = (bcf_hdr_t *) (*env)->GetLongArrayElement(env,array,1);
if(header!=NULL) bcf_hdr_destroy(header);

bcf1_t* bcf = (bcf1_t *) (*env)->GetLongArrayElement(env,array,2);
i(bcf!=NULL) bcf_destroy(bcf);
}


jint QUALIFIEDMETHOD(kseq_1read4y)(JNIEnv *env, jclass clazz,jlong hdrptr) {
bcf_hdr_t *h = (bcf_hdr_t *)hdrptr;
kstring_t htxt = {0,0,0};
if (bcf_hdr_format(h, 1, &htxt) < 0) {
     free(htxt.s);
    return -1;
}
kputc('\0', &htxt); // include the \0 byte
jobject s =  (*env)->NewStringUTF(env,htxt.s);
free(htxt.s);
return s;
}

jint QUALIFIEDMETHOD(kseq_1read4y)(JNIEnv *env, jclass clazz,jlong fpptr,jlong hdrptr,jlong bcfptr) {
{
htsFile *fp = (htsFile*)fpptr;
bcf_hdr_t *h = (bcf_hdr_t *)hdrptr;
bcf1_t* v = (bcf1_t *) bcfptr;
int rez= bcf_read(fp,h,v);
if (vcf_format1(h, v, &fp->line) != 0)
        return -1;
jobject s =  (*env)->NewStringUTF(env,fp->line.s);
return s;
}

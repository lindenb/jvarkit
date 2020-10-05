

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


typedef struct VcfShuttle {
    htsFile *in;
    bcf_hdr_t *header;
    bcf1_t* bcf;
    //
    tbx_t *tbx_idx;
    hts_idx_t *bcf_idx;
    //
    hts_itr_t *itr;
} VcfStuff ;

void VcfShuttleDestroy(VcfShuttle* vcf) {
if (vcf==NULL) return;
if (vcf->in!=NULL) hts_close(vcf->in);
if (vcf->header!=NULL) bcf_hdr_destroy(vcf->header);
if (vcf->bcf!=NULL) bcf_destroy(vcf->bcf);
if (vcf->itr!=NULL) hts_itr_destroy(vcf->itr);
if ( vcf->tbx_idx!=NULL ) tbx_destroy(vcf->tbx_idx);
if ( vcf->bcf_idx!=NULL ) hts_idx_destroy(vcf->bcf_idx);
free(vcf);
}

jlong QUALIFIEDMETHOD(kseq_1read4y)(JNIEnv *env, jclass clazz,jlong ptr,jbool require_index) {
  VcfShuttle* vcf  = (VcfShuttle*)calloc(1,sizeof(VcfShuttle));
  const char* fn = (*env)->GetStringUTFChars(env,(jstring) filename, NULL);
  vcf->in = hts_open(fn,"r");
  if(vcf->in==NULL) {
    ERROR("Cannot open input vcf %s. (%s)",fn,strerror(errno));
    VcfShuttleDestroy(vcf);    
    return 0L;    
    }
 vcf->header = bcf_hdr_read(in);
 if(vcf->header==NULL) {
    ERROR("Cannot read header for input vcf %s.",fn);
    VcfShuttleDestroy(vcf);    
    return 0L;    
    }

 if(require_index) {
    if ( vcf->in->format.format==vcf ) {
        vcf->tbx_idx = tbx_index_load(fn);
        if(vcf->tbx_idx==NULL) {
            ERROR("Cannot read index for input vcf %s.",fn);
            VcfShuttleDestroy(vcf);    
            return 0L;    
            }
        }
    else  if ( vcf->in->format.format==bcf ) {
        reader->bcf_idx = bcf_index_load(fn);
        if(vcf->bcf_idx==NULL) {
            ERROR("Cannot read index for input bcf %s.",fn);
            VcfShuttleDestroy(vcf);    
            return 0L;    
            }
        }
    else {
        ERROR("Unknown filetype for input vcf %s.",fn);
        VcfShuttleDestroy(vcf);    
        return 0L;
        }
    }

 vcf->bcf = bcf_init();
 if(vcf->bcf==NULL) {
    ERROR("Out of memory.");
    VcfShuttleDestroy(vcf);    
    return -1;    
    }

  return (jlong)vcf;
  }

jint QUALIFIEDMETHOD(kseq_1read4y)(JNIEnv *env, jclass clazz,jlong ptr) {
VcfShuttle* vcf  = (VcfShuttle*)ptr;
VcfShuttleDestroy(ptr);
return 0;
}


jlong QUALIFIEDMETHOD(kseq_1read4y)(JNIEnv *env, jclass clazz,jlong hdrptr) {
VcfShuttle* vcf  = (VcfShuttle*)ptr;
kstring_t htxt = {0,0,0};
if (bcf_hdr_format(vcf->header, 1, &htxt) < 0) {
     free(htxt.s);
    return -1;
}
kputc('\0', &htxt); // include the \0 byte
jobject s =  (*env)->NewStringUTF(env,htxt.s);
free(htxt.s);
return  (jlong)s;
}

jlong QUALIFIEDMETHOD(kseq_1read4y)(JNIEnv *env, jclass clazz,jlong ptr) {
{
VcfShuttle* vcf  = (VcfShuttle*)ptr;
int rez= bcf_read(vcf->in,vcf->header,vcf->bcf);
if (vcf_format1(vcf->header, vcf->bcf, &(vcf->in->line)) != 0)
        return -1;
jobject s =  (*env)->NewStringUTF(env,vcf->in->line.s);
return (jlong)s;
}

jint QUALIFIEDMETHOD(kseq_1read4y)(JNIEnv *env, jclass clazz,jlong ptr) {
VcfShuttle* vcf  = (VcfShuttle*)ptr;
if ( vcf->itr!=NULL )
    {
        hts_itr_destroy(vcf->itr);
        vcf->itr = NULL;
    }
if ( vcf->tbx_idx !=NULL)
    {
        int tid = tbx_name2id(vcf->tbx_idx, seq);
        if ( tid==-1 ) return -1;    // the sequence not present in this file
        vcf->itr = tbx_itr_queryi(vcf->tbx_idx,tid,start,end+1);
    }
else
    {
        int tid = bcf_hdr_name2id(vcf->header, seq);
        if ( tid==-1 ) return -1;    // the sequence not present in this file
        vcf->itr = bcf_itr_queryi(vcf->bcf_idx,tid,start,end+1);
    }
}




 int ret = 0;
 if ( reader->tbx_idx )
        {
            if ( (ret=tbx_itr_next(vcf->in, vcf->tbx_idx, vcf->itr, &files->tmps)) < 0 ) break;  // no more lines
            ret = vcf_parse1(&files->tmps, vcf->header, reader->buffer[reader->nbuffer+1]);
            if ( ret<0 ) { files->errnum = vcf_parse_error; break; }
        }
else
        {
            ret = bcf_itr_next(reader->file, reader->itr, reader->buffer[reader->nbuffer+1]);
            if ( ret < -1 ) files->errnum = bcf_read_error;
            if ( ret < 0 ) break; // no more lines or an error
            bcf_subset_format(reader->header,reader->buffer[reader->nbuffer+1]);
        }
if ( files->streaming )
        {
            if ( reader->file->format.format==vcf )
            {
                if ( (ret=hts_getline(reader->file, KS_SEP_LINE, &files->tmps)) < 0 ) break;   // no more lines
                ret = vcf_parse1(&files->tmps, reader->header, reader->buffer[reader->nbuffer+1]);
                if ( ret<0 ) { files->errnum = vcf_parse_error; break; }
            }
            else if ( reader->file->format.format==bcf )
            {
                ret = bcf_read1(reader->file, reader->header, reader->buffer[reader->nbuffer+1]);
                if ( ret < -1 ) files->errnum = bcf_read_error;
                if ( ret < 0 ) break; // no more lines or an error
            }
            else
            {
                hts_log_error("Fixme: not ready for this");
                exit(1);
            }
        }


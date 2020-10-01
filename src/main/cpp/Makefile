.SHELL=/bin/bash
HTSLIB_VERSION?=1.11
CC?=gcc
CFLAGS=-O3 -Wall  -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE
ARCH=$(shell arch)

ifeq (${JAVA_HOME},)
$(error $${JAVA_HOME} is not defined)
endif

## find path where to find include files
JDK_JNI_INCLUDES?=$(addprefix -I,$(sort $(dir $(shell find ${JAVA_HOME}/include -type f -name "*.h"))))


ifeq (${JDK_JNI_INCLUDES},)
$(error Cannot find C header files under $${JAVA_HOME})
endif

## see https://github.com/lindenb/jbwa/pull/5
ifeq (${OSNAME},Darwin)
native.extension=jnilib
else
native.extension=so
endif



libhtsjni.$(ARCH).${native.extension} : com_github_lindenb_jvarkit_htslib_HtsLib.o htslib-$(HTSLIB_VERSION)/libhts.a
	$(CC) -dynamiclib -shared -o $@ $<  -L htslib-$(HTSLIB_VERSION) -lhts -lm -lz -lpthread

com_github_lindenb_jvarkit_htslib_HtsLib.o : com_github_lindenb_jvarkit_htslib_HtsLib.c com_github_lindenb_jvarkit_htslib_HtsLib.h htslib-$(HTSLIB_VERSION)/libhts.a
	$(CC) -c $(CFLAGS) -o $@ -fPIC  $(JDK_JNI_INCLUDES)  -I htslib-$(HTSLIB_VERSION) $<

com_github_lindenb_jvarkit_htslib_HtsLib.h : ../java/com/github/lindenb/jvarkit/htslib/HtsLib.java
	javah -o $@ -jni -force -classpath  ../java "com.github.lindenb.jvarkit.htslib.HtsLib"

htslib-$(HTSLIB_VERSION)/libhts.a :  htslib-$(HTSLIB_VERSION)/Makefile
	cd $(dir $@) && make


htslib-$(HTSLIB_VERSION)/Makefile:
	rm -rf "$(dir $@)"
	wget -O "${HTSLIB_VERSION}.tar.bz2" "https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2" || curl -o "${HTSLIB_VERSION}.tar.bz2" "https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2"
	tar xfj "${HTSLIB_VERSION}.tar.bz2"
	rm -vf "${HTSLIB_VERSION}.tar.bz2"
	touch -c $@

clean:
	rm -rf htslib-$(HTSLIB_VERSION) lib

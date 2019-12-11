ARCH=$(shell uname -m)
HTSLIB_VERSION=1.10
CC=gcc

## find path where to find include files
JDK_JNI_INCLUDES?=$(addprefix -I,$(sort $(dir $(shell find ${JAVA_HOME}/include -type f -name "*.h"))))


ifeq (${JDK_JNI_INCLUDES},)
$(error Cannot find C header files under $${JAVA_HOME})
endif


lib/jvarkit.$(ARCH).a : htslib/libhts.a
	mkdir -p "$(dir $@)"
	# https://stackoverflow.com/questions/915128
	$(CC) -shared -o $@ -Wl,--whole-archive $(filter %.a,$^) -Wl,--no-whole-archive

htslib/libhts.a :  htslib/Makefile
	cd $(dir $@) && make
htslib/Makefile:
	rm -rf "$(dir $@)"
	wget -O "${HTSLIB_VERSION}.tar.gz" "https://github.com/samtools/htslib/archive/${HTSLIB_VERSION}.tar.gz"
	tar xfz "${HTSLIB_VERSION}.tar.gz"
	mv -v "htslib-${HTSLIB_VERSION}" "$(dir $@)"
	rm -vf "${HTSLIB_VERSION}.tar.gz"
	touch -c $@

clean:
	rm -rf htslib lib

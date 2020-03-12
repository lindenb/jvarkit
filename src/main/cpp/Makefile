

htslib/libhts.a :  htslib/Makefile
	cd $(dir $@) && make
htslib/Makefile:
	rm -rf "$(dir $@)"
	wget -O "${HTSLIB_VERSION}.tar.bz2" "https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2"
	tar xfj "${HTSLIB_VERSION}.tar.bz2"
	mv -v "htslib-${HTSLIB_VERSION}" "$(dir $@)"
	rm -vf "${HTSLIB_VERSION}.tar.bz2"
	touch -c $@

clean:
	rm -rf htslib lib

FROM ubuntu:22.04

LABEL about.summary = "Jvarkit: java utilities for bioinformatics https://github.com/lindenb/jvarkit. The jar archive for jvarkit is located in  /opt/jvarkit/dist/jvarkit.jar. This container also includes R, datamash, bwa, samtools, htslib, bcftools, bedtools and /opt/picard/picard.jar"
LABEL about.home = "https://github.com/lindenb/jvarkit"
LABEL about.documentation = "https://jvarkit.readthedocs.io/en/latest/"
LABEL about.tag = "Bioinformatics,NGS,Genomics,BAM,VCF"

MAINTAINER Pierre Lindenbaum PhD Institut-du-Thorax Nantes/France

ARG JVARKIT_VERSION=1b2aedf24
ARG HTS_VERSION=1.17
ARG BEDTOOLS_VERSION=2.30.0
ARG BWA_VERSION=139f68fc4c37478137
ARG PICARD_VERSION=2.27.5
ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get -y update && \
	apt-get -y install git wget openjdk-17-jdk make r-base bzip2 libcurl4-openssl-dev libncurses5-dev python3 libz-dev libbz2-dev liblzma-dev datamash && \
	mkdir -p /opt/picard && \
	wget -O /opt/picard/picard.jar "https://github.com/broadinstitute/picard/releases/download/${PICARD_VERSION}/picard.jar" && \
	cd /opt/ && \
	git clone -b dev "https://github.com/lindenb/jvarkit.git" jvarkit.tmp && \
	(cd jvarkit.tmp && git reset --hard ${JVARKIT_VERSION} && ./gradlew  -Djvarkit.disable.test=true jvarkit) && \
	mv "jvarkit.tmp" "/opt/jvarkit" && \
    (cd /opt && wget -O bedtools-${BEDTOOLS_VERSION}.tar.gz "https://github.com/arq5x/bedtools2/releases/download/v${BEDTOOLS_VERSION}/bedtools-${BEDTOOLS_VERSION}.tar.gz" && tar xvfz bedtools-${BEDTOOLS_VERSION}.tar.gz && cd bedtools2 && sed -i 's/python scripts/python3 scripts/' Makefile && make) && \
    (cd /opt/ && wget -O htslib-${HTS_VERSION}.tar.bz2 "https://github.com/samtools/htslib/releases/download/${HTS_VERSION}/htslib-${HTS_VERSION}.tar.bz2" && tar xvfj htslib-${HTS_VERSION}.tar.bz2 && cd htslib-${HTS_VERSION} && make) && \
    (cd /opt/ && wget -O bcftools-${HTS_VERSION}.tar.bz2 "https://github.com/samtools/bcftools/releases/download/${HTS_VERSION}/bcftools-${HTS_VERSION}.tar.bz2" && tar xvfj bcftools-${HTS_VERSION}.tar.bz2 && cd bcftools-${HTS_VERSION} && make HTSDIR=/opt/htslib-${HTS_VERSION}) && \
     (cd /opt/ && wget -O samtools-${HTS_VERSION}.tar.bz2 "https://github.com/samtools/samtools/releases/download/${HTS_VERSION}/samtools-${HTS_VERSION}.tar.bz2" && tar xvfj samtools-${HTS_VERSION}.tar.bz2 && cd samtools-${HTS_VERSION} && make HTSDIR=/opt/htslib-${HTS_VERSION}) && \
    (cd /opt && git clone "https://github.com/lh3/bwa" && cd bwa && git reset --hard "${BWA_VERSION}" && make) && \
    rm -f /opt/*.tar.bz2 /opt/*.tar.gz

RUN apt-get clean

ENV JVARKIT_DIST="/opt/jvarkit/dist"
ENV JVARKIT_JAR="/opt/jvarkit/dist/jvarkit.jar"
ENV PICARD_JAR="/opt/picard/picard.jar"
ENV BCFTOOLS_PLUGINS="/opt/bcftools-${HTS_VERSION}/plugins"
ENV LD_LIBRARY_PATH="/opt/htslib-${HTS_VERSION}:${LD_LIBRARY_PATH}"
ENV PATH="/opt/bwa:/opt/bedtools2/bin:/opt/htslib-${HTS_VERSION}:/opt/bcftools-${HTS_VERSION}:/opt/samtools-${HTS_VERSION}:${PATH}"



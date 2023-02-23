FROM ubuntu:22.04

LABEL description="Jvarkit: java utilities for bioinformatics https://github.com/lindenb/jvarkit. The jar archive for jvarkit is located in  /opt/jvarkit/dist/jvarkit.jar. This container also includes bwa, samtools, bcftools, bedtools and /opt/picard/picard.jar"


MAINTAINER Pierre Lindenbaum PhD Institut-du-Thorax Nantes/France

ARG HTS_VERSION=1.17
ARG BEDTOOLS_VERSION=2.30.0
ARG BWA_VERSION=0.7.17

RUN apt-get -y update && \
	apt-get -y install git wget openjdk-11-jdk make && \
	mkdir -p /opt/picard && \
	wget -O /opt/picard/picard.jar "https://github.com/broadinstitute/picard/releases/download/2.27.5/picard.jar" && \
	cd /opt/ && \
	git clone -b dev "https://github.com/lindenb/jvarkit.git" jvarkit.tmp && \
	(cd jvarkit.tmp && git reset --hard c937fde && ./gradlew  -Djvarkit.disable.test=true jvarkit) && \
	mv "jvarkit.tmp" "/opt/jvarkit" && \
    (cd /opt/ && wget -O htslib-${HTS_VERSION}.tar.bz2 "https://github.com/samtools/htslib/releases/download/${HTS_VERSION}/htslib-${HTS_VERSION}.tar.bz2" && tar xvfj htslib-${HTS_VERSION}.tar.bz2 && cd htslib-${HTS_VERSION} && make) && \
    (cd /opt/ && wget -O bcftools-${HTS_VERSION}.tar.bz2 "https://github.com/samtools/bcftools/releases/download/${HTS_VERSION}/bcftools-${HTS_VERSION}.tar.bz2" && tar xvfj bcftools-${HTS_VERSION}.tar.bz2 && cd bcftools-${HTS_VERSION} && make HTSDIR=/opt/htslib-${HTS_VERSION}) && \
     (cd /opt/ && wget -O samtools-${HTS_VERSION}.tar.bz2 "https://github.com/samtools/samtools/releases/download/${HTS_VERSION}/samtools-${HTS_VERSION}.tar.bz2" && tar xvfj samtools-${HTS_VERSION}.tar.bz2 && cd samtools-${HTS_VERSION} && make HTSDIR=/opt/htslib-${HTS_VERSION}) && \
    (cd /opt && wget -O bedtools-${BEDTOOLS_VERSION}.tar.gz "https://github.com/arq5x/bedtools2/releases/download/v${BEDTOOLS_VERSION}/bedtools-${BEDTOOLS_VERSION}.tar.gz" && tar xvfz bedtools-${BEDTOOLS_VERSION}.tar.gz && cd bedtools2 && make) && \
    (cd /opt && wget -O bwa-${BWA_VERSION}.tar.bz2 "https://github.com/lh3/bwa/releases/download/v${BWA_VERSION}/bwa-${BWA_VERSION}.tar.bz2" && tar xvfj bwa-${BWA_VERSION}.tar.bz2 && cd  bwa-${BWA_VERSION} && make) && \
    rm -f /opt/*.tar.bz2 /opt/*.tar.gz && \
    mkdir -p /etc/profile.d && \
    echo -e 'export JVARKIT_DIST=/opt/jvarkit/dist\nexport PICARD_JAR=/opt/picard/picard.jar\nexport BCFTOOLS_PLUGINS=/opt/bcftools-${HTS_VERSION}/plugins\nexport PATH=/opt/bwa-${BWA_VERSION}:/opt/bedtools2/bin:/opt/htslib-${HTS_VERSION}:/opt/bcftools-${HTS_VERSION}:/opt/samtools-${HTS_VERSION}:$${PATH}' > /etc/profile.d/jvarkit.sh


RUN apt-get clean



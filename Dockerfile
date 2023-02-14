FROM ubuntu:22.04

LABEL description="Jvarkit: java utilities for bioinformatics https://github.com/lindenb/jvarkit. The jar archive for jvarkit is located in  /opt/jvarkit/dist/jvarkit.jar. This container also includes bwa, samtools, bcftools, bedtools and /opt/picard/picard.jar"


MAINTAINER Pierre Lindenbaum PhD Institut-du-Thorax Nantes/France

RUN apt-get -y update && \
	apt-get -y install git wget openjdk-11-jdk make samtools bcftools bedtools bwa && \
	mkdir -p /opt/picard && \
	wget -O /opt/picard/picard.jar "https://github.com/broadinstitute/picard/releases/download/2.27.5/picard.jar" && \
	cd /opt/ && \
	git clone -b dev "https://github.com/lindenb/jvarkit.git" jvarkit.tmp && \
	(cd jvarkit.tmp && git reset --hard 1e09f06d4c05e5a148 && ./gradlew  -Djvarkit.disable.test=true jvarkit) && \
	mv "jvarkit.tmp" "/opt/jvarkit" 

ENV JVARKIT_DIST /opt/jvarkit/dist
RUN export JVARKIT_DIST


RUN apt-get clean



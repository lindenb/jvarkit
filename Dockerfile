FROM ubuntu:22.04

LABEL Jvarkit: java utilities for bioinformatics https://github.com/lindenb/jvarkit


MAINTAINER Pierre Lindenbaum PhD Institut-du-Thorax Nantes/France

RUN apt-get -y update && \
	apt-get -y install git wget openjdk-11-jdk make && \
	mkdir -p /opt/ && \
	cd /opt/ && \
	git clone -b dev "https://github.com/lindenb/jvarkit.git" jvarkit.tmp && \
	(cd jvarkit.tmp && git reset --hard 9a2f3d5ac28e22802 && ./gradlew  -Djvarkit.disable.test=true jvarkit) && \
	mv "jvarkit.tmp" "/opt/jvarkit" 

ENV JVARKIT_DIST /opt/jvarkit/dist
RUN export JVARKIT_DIST


RUN apt-get clean



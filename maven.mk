#
# libraries to be downloaded by maven:
#
lib.dir?=lib

avro.tools.version = 1.7.7
avro.libs = $(lib.dir)/org/apache/avro/avro-tools/${avro.tools.version}/avro-tools-${avro.tools.version}.jar

commons.loggging.jars = \
	$(lib.dir)/commons-logging/commons-logging/1.2/commons-logging-1.2.jar

httpclient.libs  = \
	${commons.loggging.jars} \
	$(lib.dir)/org/apache/httpcomponents/httpcore/4.4.1/httpcore-4.4.1.jar \
	$(lib.dir)/org/apache/httpcomponents/httpclient/4.5/httpclient-4.5.jar \
	$(lib.dir)/commons-codec/commons-codec/1.10/commons-codec-1.10.jar
	

common.math3.libs  =  \
	$(lib.dir)/org/apache/commons/commons-math3/3.5/commons-math3-3.5.jar

apache.commons.cli.jars  = \
	$(lib.dir)/commons-cli/commons-cli/1.3.1/commons-cli-1.3.1.jar

commons.beanutils.jars  =  \
	${commons.loggging.jars} \
	$(lib.dir)/commons-collections/commons-collections/3.2.1/commons-collections-3.2.1.jar \
	$(lib.dir)/commons-beanutils/commons-beanutils/1.9.2/commons-beanutils-1.9.2.jar

commons.validator.jars  = \
	${commons.beanutils.jars} \
	$(lib.dir)/commons-validator/commons-validator/1.4.1/commons-validator-1.4.1.jar \
	$(lib.dir)/commons-digester/commons-digester/1.8.1/commons-digester-1.8.1.jar 

slf4j.jars = \
	$(lib.dir)/org/slf4j/slf4j-api/1.7.13/slf4j-api-1.7.13.jar \
	$(lib.dir)/org/slf4j/slf4j-simple/1.7.13/slf4j-simple-1.7.13.jar

all_maven_jars = $(sort  ${slf4j.jars} ${httpclient.libs} ${avro.libs} ${common.math3.libs} ${apache.commons.cli.jars} ${commons.validator.jars})

${all_maven_jars}  : 
	mkdir -p $(dir $@) && curl -Lk ${curl.proxy} -o "$@" "http://central.maven.org/maven2/$(patsubst ${lib.dir}/%,%,$@)"


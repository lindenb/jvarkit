#
# libraries to be downloaded by maven:
#
lib.dir?=lib

avro.tools.version = 1.7.7
avro.libs = $(lib.dir)/org/apache/avro/avro-tools/${avro.tools.version}/avro-tools-${avro.tools.version}.jar


lib.dir?=lib

htsjdk.version=2.9.1
htsjdk.jars  =  \
	$(lib.dir)/com/github/samtools/htsjdk/${htsjdk.version}/htsjdk-${htsjdk.version}.jar \
	$(lib.dir)/commons-logging/commons-logging/1.1.1/commons-logging-1.1.1.jar \
	$(lib.dir)/gov/nih/nlm/ncbi/ngs-java/1.2.4/ngs-java-1.2.4.jar \
	$(lib.dir)/org/apache/commons/commons-compress/1.4.1/commons-compress-1.4.1.jar \
	$(lib.dir)/org/apache/commons/commons-jexl/2.1.1/commons-jexl-2.1.1.jar \
	$(lib.dir)/org/tukaani/xz/1.5/xz-1.5.jar \
	$(lib.dir)/org/xerial/snappy/snappy-java/1.0.3-rc3/snappy-java-1.0.3-rc3.jar

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

derby.jars  =  \
	$(lib.dir)/org/apache/derby/derby/10.12.1.1/derby-10.12.1.1.jar \
	$(lib.dir)/org/apache/derby/derbyclient/10.12.1.1/derbyclient-10.12.1.1.jar

derby-tools.jar = \
	$(lib.dir)/org/apache/derby/derbytools/10.12.1.1/derbytools-10.12.1.1.jar \

jetty.jars  =  \
	$(lib.dir)/javax/servlet/javax.servlet-api/4.0.0-b01/javax.servlet-api-4.0.0-b01.jar \
	$(lib.dir)/org/eclipse/jetty/jetty-webapp/9.3.7.v20160115/jetty-webapp-9.3.7.v20160115.jar \
	$(lib.dir)/org/eclipse/jetty/jetty-http/9.3.7.v20160115/jetty-http-9.3.7.v20160115.jar \
	$(lib.dir)/org/eclipse/jetty/jetty-server/9.3.7.v20160115/jetty-server-9.3.7.v20160115.jar \
	$(lib.dir)/org/eclipse/jetty/jetty-io/9.3.7.v20160115/jetty-io-9.3.7.v20160115.jar \
	$(lib.dir)/org/eclipse/jetty/jetty-security/9.3.7.v20160115/jetty-security-9.3.7.v20160115.jar \
	$(lib.dir)/org/eclipse/jetty/jetty-servlet/9.3.7.v20160115/jetty-servlet-9.3.7.v20160115.jar \
	$(lib.dir)/org/eclipse/jetty/jetty-util/9.3.7.v20160115/jetty-util-9.3.7.v20160115.jar \
	$(lib.dir)/org/eclipse/jetty/jetty-xml/9.3.7.v20160115/jetty-xml-9.3.7.v20160115.jar

spring-beans.jars  =  \
	$(lib.dir)/org/springframework/spring-core/4.2.4.RELEASE/spring-core-4.2.4.RELEASE.jar \
	$(lib.dir)/org/springframework/spring-beans/4.2.4.RELEASE/spring-beans-4.2.4.RELEASE.jar \
	$(lib.dir)/commons-logging/commons-logging/1.2/commons-logging-1.2.jar

web.frameworks.jar  =  \
	$(lib.dir)/org/webjars/bootstrap/3.3.5/bootstrap-3.3.5.jar \
	$(lib.dir)/org/webjars/jquery-ui/1.11.4/jquery-ui-1.11.4.jar \
	$(lib.dir)/org/webjars/jquery/1.11.1/jquery-1.11.1.jar


gson.jar = \
	$(lib.dir)/com/google/code/gson/gson/2.6.2/gson-2.6.2.jar

velocity.jars  =  \
	$(lib.dir)/commons-collections/commons-collections/3.2.1/commons-collections-3.2.1.jar \
	$(lib.dir)/commons-lang/commons-lang/2.4/commons-lang-2.4.jar \
	$(lib.dir)/org/apache/velocity/velocity/1.7/velocity-1.7.jar

jcommander.jar= \
	$(lib.dir)/com/beust/jcommander/1.64/jcommander-1.64.jar


javacc.jar=\
	$(lib.dir)/net/java/dev/javacc/javacc/7.0.2/javacc-7.0.2.jar

berkeleydb.jar=$(lib.dir)/com/sleepycat/je/7.3.7/je-7.3.7.jar
${berkeleydb.jar}:
	mkdir -p $(dir $@) && wget -O "$@" "http://download.oracle.com/maven/$(patsubst ${lib.dir}/%,%,$@)"

drools.jar  =  \
	$(lib.dir)/com/google/protobuf/protobuf-java/3.3.1/protobuf-java-3.3.1.jar \
	$(lib.dir)/com/thoughtworks/xstream/xstream/1.4.10-java7/xstream-1.4.10-java7.jar \
	$(lib.dir)/commons-codec/commons-codec/1.10/commons-codec-1.10.jar \
	$(lib.dir)/javax/activation/activation/1.1.1/activation-1.1.1.jar \
	$(lib.dir)/jmock/jmock/1.0.0/jmock-1.0.0.jar \
	$(lib.dir)/org/antlr/antlr-runtime/3.5.2/antlr-runtime-3.5.2.jar \
	$(lib.dir)/org/drools/drools-compiler/7.0.0.Final/drools-compiler-7.0.0.Final.jar \
	$(lib.dir)/org/drools/drools-core/7.0.0.Final/drools-core-7.0.0.Final.jar \
	$(lib.dir)/org/eclipse/jdt/core/compiler/ecj/4.6.1/ecj-4.6.1.jar \
	$(lib.dir)/org/hamcrest/hamcrest-core/1.3/hamcrest-core-1.3.jar \
	$(lib.dir)/org/kie/kie-api/7.1.0.Beta2/kie-api-7.1.0.Beta2.jar \
	$(lib.dir)/org/kie/kie-internal/7.1.0.Beta2/kie-internal-7.1.0.Beta2.jar \
	$(lib.dir)/org/mvel/mvel2/2.3.1.Final/mvel2-2.3.1.Final.jar \
	$(lib.dir)/org/slf4j/slf4j-api/1.8.0-alpha2/slf4j-api-1.8.0-alpha2.jar \
	$(lib.dir)/xmlpull/xmlpull/1.1.2.1/xmlpull-1.1.2.1.jar \
	$(lib.dir)/xpp3/xpp3_min/1.1.3.4.O/xpp3_min-1.1.3.4.O.jar

	
all_maven_jars = $(sort ${drools.jar} ${javacc.jar} ${jcommander.jar} ${velocity.jars} ${htsjdk.jars} ${web.frameworks.jar} ${spring-beans.jars} ${jetty.jars} ${derby.jars} ${slf4j.jars} ${httpclient.libs} ${avro.libs} ${common.math3.libs} ${apache.commons.cli.jars} ${commons.validator.jars} ${gson.jar} ${derby-tools.jar} )


${all_maven_jars}  : 
	mkdir -p $(dir $@) && curl -Lk ${curl.proxy} -o "$@" "http://central.maven.org/maven2/$(patsubst ${lib.dir}/%,%,$@)"


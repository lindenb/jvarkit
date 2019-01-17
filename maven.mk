#
# libraries to be downloaded by maven:
#
lib.dir?=lib

avro.tools.version = 1.7.7
avro.libs = $(lib.dir)/org/apache/avro/avro-tools/${avro.tools.version}/avro-tools-${avro.tools.version}.jar


lib.dir?=lib

mysql.jar = \
	 $(lib.dir)/mysql/mysql-connector-java/5.1.47/mysql-connector-java-5.1.47.jar

htsjdk.version=2.18.2
htsjdk.jars  =  \
	$(lib.dir)/com/github/samtools/htsjdk/2.18.2/htsjdk-2.18.2.jar \
	$(lib.dir)/commons-logging/commons-logging/1.1.1/commons-logging-1.1.1.jar \
	$(lib.dir)/gov/nih/nlm/ncbi/ngs-java/2.9.0/ngs-java-2.9.0.jar \
	$(lib.dir)/org/apache/commons/commons-compress/1.4.1/commons-compress-1.4.1.jar \
	$(lib.dir)/org/apache/commons/commons-jexl/2.1.1/commons-jexl-2.1.1.jar \
	$(lib.dir)/org/tukaani/xz/1.5/xz-1.5.jar \
	$(lib.dir)/org/xerial/snappy/snappy-java/1.1.4/snappy-java-1.1.4.jar



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
	$(lib.dir)/org/springframework/spring-context/4.2.4.RELEASE/spring-context-4.2.4.RELEASE.jar \
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


testng.jars = \
	$(lib.dir)/org/testng/testng/6.14.3/testng-6.14.3.jar \
	${jcommander.jar}

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

spring.batch.jars = \
	$(lib.dir)/com/ibm/jbatch/com.ibm.jbatch-tck-spi/1.0/com.ibm.jbatch-tck-spi-1.0.jar \
	$(lib.dir)/com/thoughtworks/xstream/xstream/1.4.7/xstream-1.4.7.jar \
	$(lib.dir)/commons-logging/commons-logging/1.1.3/commons-logging-1.1.3.jar \
	$(lib.dir)/javax/batch/javax.batch-api/1.0/javax.batch-api-1.0.jar \
	$(lib.dir)/org/codehaus/jettison/jettison/1.2/jettison-1.2.jar \
	$(lib.dir)/org/springframework/batch/spring-batch-core/3.0.8.RELEASE/spring-batch-core-3.0.8.RELEASE.jar \
	$(lib.dir)/org/springframework/batch/spring-batch-infrastructure/3.0.8.RELEASE/spring-batch-infrastructure-3.0.8.RELEASE.jar \
	$(lib.dir)/org/springframework/retry/spring-retry/1.2.1.RELEASE/spring-retry-1.2.1.RELEASE.jar \
	$(lib.dir)/org/springframework/spring-aop/5.0.0.RELEASE/spring-aop-5.0.0.RELEASE.jar \
	$(lib.dir)/org/springframework/spring-beans/5.0.0.RELEASE/spring-beans-5.0.0.RELEASE.jar \
	$(lib.dir)/org/springframework/spring-context/5.0.0.RELEASE/spring-context-5.0.0.RELEASE.jar \
	$(lib.dir)/org/springframework/spring-core/5.0.0.RELEASE/spring-core-5.0.0.RELEASE.jar \
	$(lib.dir)/org/springframework/spring-expression/5.0.0.RELEASE/spring-expression-5.0.0.RELEASE.jar \
	$(lib.dir)/org/springframework/spring-tx/4.0.5.RELEASE/spring-tx-4.0.5.RELEASE.jar

jaxb.jars = \
	$(lib.dir)/javax/xml/bind/jaxb-api/2.2.11/jaxb-api-2.2.11.jar \
	$(lib.dir)/javax/activation/javax.activation-api/1.2.0/javax.activation-api-1.2.0.jar \
	$(lib.dir)/org/glassfish/jaxb/jaxb-core/2.2.11/jaxb-core-2.2.11.jar \
	$(lib.dir)/org/glassfish/jaxb/jaxb-runtime/2.2.11/jaxb-runtime-2.2.11.jar\
	$(lib.dir)/com/sun/istack/istack-commons-runtime/3.0.8/istack-commons-runtime-3.0.8.jar
	

	
	
JAXB_VERSION=2.2.11

xjc_jars= \
	$(lib.dir)/com/beust/jcommander/1.72/jcommander-1.72.jar \
	$(lib.dir)/com/sun/istack/istack-commons-runtime/3.0.8/istack-commons-runtime-3.0.8.jar \
	$(lib.dir)/com/sun/istack/istack-commons-tools/3.0.8/istack-commons-tools-3.0.8.jar \
	$(lib.dir)/com/sun/xml/bind/external/rngom/2.2.11/rngom-2.2.11.jar \
	$(lib.dir)/com/sun/xml/dtd-parser/dtd-parser/1.4.1/dtd-parser-1.4.1.jar \
	$(lib.dir)/com/sun/xsom/xsom/20140925/xsom-20140925.jar \
	$(lib.dir)/jakarta/activation/jakarta.activation-api/1.2.1/jakarta.activation-api-1.2.1.jar \
	$(lib.dir)/javax/activation/javax.activation-api/1.2.0/javax.activation-api-1.2.0.jar \
	$(lib.dir)/javax/xml/bind/jaxb-api/2.4.0-b180830.0359/jaxb-api-2.4.0-b180830.0359.jar \
	$(lib.dir)/junit/junit/4.13-beta-1/junit-4.13-beta-1.jar \
	$(lib.dir)/org/apache/ant/ant-launcher/1.10.5/ant-launcher-1.10.5.jar \
	$(lib.dir)/org/apache/ant/ant/1.10.5/ant-1.10.5.jar \
	$(lib.dir)/org/glassfish/jaxb/codemodel/2.2.11/codemodel-2.2.11.jar \
	$(lib.dir)/org/glassfish/jaxb/jaxb-core/2.3.0.1/jaxb-core-2.3.0.1.jar \
	$(lib.dir)/org/glassfish/jaxb/jaxb-xjc/2.2.11/jaxb-xjc-2.2.11.jar \
	$(lib.dir)/org/glassfish/jaxb/txw2/2.3.0.1/txw2-2.3.0.1.jar \
	$(lib.dir)/org/hamcrest/hamcrest-core/2.1/hamcrest-core-2.1.jar \
	$(lib.dir)/org/hamcrest/hamcrest/2.1/hamcrest-2.1.jar \
	$(lib.dir)/org/testng/testng/7.0.0-beta3/testng-7.0.0-beta3.jar \
	$(lib.dir)/relaxngDatatype/relaxngDatatype/20020414/relaxngDatatype-20020414.jar \
	$(lib.dir)/com/sun/xml/bind/jaxb-impl/2.2.11/jaxb-impl-2.2.11.jar 


all_maven_jars = $(sort ${mysql.jar} ${testng.jars} ${drools.jar} ${javacc.jar} ${jcommander.jar} ${velocity.jars} ${htsjdk.jars} ${web.frameworks.jar} ${spring-beans.jars} ${jetty.jars} ${derby.jars} ${slf4j.jars} ${httpclient.libs} ${avro.libs} ${common.math3.libs} ${apache.commons.cli.jars} ${commons.validator.jars} ${gson.jar} ${derby-tools.jar} ${spring.batch.jars} ${jaxb.jars} ${xjc_jars})


${all_maven_jars}  : 
	mkdir -p $(dir $@) && wget -O "$(addsuffix .tmp,$@)" "http://central.maven.org/maven2/$(patsubst ${lib.dir}/%,%,$@)" && mv -v "$(addsuffix .tmp,$@)" "$@"


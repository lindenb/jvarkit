#
# libraries to be downloaded by maven:
#
avro.tools.version = 1.7.7
avro.libs = lib/org/apache/avro/avro-tools/${avro.tools.version}/avro-tools-${avro.tools.version}.jar

${avro.libs} :
	mkdir -p $(dir $@) && curl  -Lk ${curl.proxy} -o $@ "http://central.maven.org/maven2/$(patsubst lib/%,%,$@)"


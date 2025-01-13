/*
The MIT License (MIT)

Copyright (c) 2025 Pierre Lindenbaum

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/
package com.github.lindenb.jvarkit.tools.vcfspring;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.variant.VariantAnnotator;
import com.github.lindenb.jvarkit.variant.VcfSpringBeanVariantAnnotator;
import com.github.lindenb.jvarkit.variant.vcf.AbstractOnePassVcfAnnotator;
/**

BEGIN_DOC

# Motivation

Use the spring framework ( https://docs.spring.io/spring-framework/docs/4.2.x/spring-framework-reference/html/xsd-configuration.html )  to load VCF filters implementing `com.github.lindenb.jvarkit.variant.VariantAnnotator`


# XML Configuration example

```xml
<?xml version="1.0" encoding="UTF-8"?>
<beans xmlns="http://www.springframework.org/schema/beans"
    xmlns:context="http://www.springframework.org/schema/context"
    xmlns:util="http://www.springframework.org/schema/util"
    xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
    xsi:schemaLocation="
        http://www.springframework.org/schema/beans
        http://www.springframework.org/schema/beans/spring-beans.xsd
        http://www.springframework.org/schema/util
        https://www.springframework.org/schema/util/spring-util.xsd
        "
	>
    <util:list id="main">
    	<bean class=" com.github.lindenb.jvarkit.util.vcf.JexlVariantPredicate" >
    		<constructor-arg index="0" type="String" value="!vc.isSNP()"/>
    		 <property name="softFiltering" value="true"/>
    		 <property name="filterName" value="HELLO"/>
    	</bean>
    	<bean class=" com.github.lindenb.jvarkit.util.vcf.JexlVariantPredicate" >
    		<constructor-arg index="0" type="String" value="vc.homVarCount &gt; 1"/>
    		 <property name="softFiltering" value="true"/>
    		 <property name="filterName" value="HELLO2"/>
    	</bean>
    </util:list>
</beans>
```

## Example:

```
$ java -jar dist/jvarkit.jar vcfspringfilter -c config.xml  src/test/resources/spring-variant-annotators.01.xml
##fileformat=VCFv4.2
(...)
##FILTER=<ID=HELLO,Description="Filtered with JEXL expression(s) ( A Java EXpression Language (JEXL) expressions to filter the variants from a VCF. Empty string will accept all variants. Expression returning a TRUE will accept the variant. See https://gatk.broadinstitute.org/hc/en-us/articles/360035891011  ) : !vc.isSNP()">
##FILTER=<ID=HELLO2,Description="Filtered with JEXL expression(s) ( A Java EXpression Language (JEXL) expressions to filter the variants from a VCF. Empty string will accept all variants. Expression returning a TRUE will accept the variant. See https://gatk.broadinstitute.org/hc/en-us/articles/360035891011  ) : vc.homVarCount > 1">
(...)
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4	S5
(...)
RF03	2573	.	A	G	17.83	HELLO	AC=4;AN=10;BQB=0.974597;DP=9;DP4=0,5,0,4;HOB=0.48;ICB=0.117361;MQ=60;MQ0F=0;MQB=0.974597;RPB=0.487298;SGB=0.473945;VDB=0.0516381	GT:PL	0/0:0,3,17	1/1:31,6,0	1/1:31,6,0	0/0:0,3,17	0/0:0,9,42
RF04	887	.	A	G	5.31	HELLO;HELLO2	AC=1;AN=10;BQB=1;DP=48;DP4=16,28,3,1;HOB=0.02;ICB=0.0439024;MQ=60;MQ0F=0;MQB=1;MQSB=1;RPB=0.90467;SGB=3.91248;VDB=0.811811	GT:PL	0/1:40,0,28	0/0:0,24,98	0/0:0,24,98	0/0:0,33,120	0/0:0,42,134
RF08	992	.	G	C	70	HELLO	AC=4;AN=10;BQB=1;DP=33;DP4=0,21,0,12;HOB=0.48;ICB=0.117361;MQ=60;MQ0F=0;MQB=1;RPB=0.73431;SGB=7.42075;VDB=0.750182	GT:PL	0/0:0,21,66	1/1:62,18,0	1/1:62,18,0	0/0:0,15,57	0/0:0,27,72
(...)
```


END_DOC

**/
@Program(
	name="vcfspringfilter",
	description="Uses the java spring Framework to build complex vcf filters",
	keywords={"vcf","java","spring","framework"},
	jvarkit_amalgamion = true,
	creationDate = "20230526",
	modificationDate = "20230526",
	menu="VCF Manipulation"
	)
public class VcfSpringFilter extends AbstractOnePassVcfAnnotator {
	private static final Logger LOG=Logger.build(VcfSpringFilter.class).make();

	@Parameter(names={"-c","--config"},description="Spring XML configuration file (  https://docs.spring.io/spring-framework/docs/4.2.x/spring-framework-reference/html/xsd-configuration.html ) .", required = true)
	private List<Path> springCongigFiles = new ArrayList<>();
	@Parameter(names={"-m","--main"},description="Main bean name")
	private String mainBeanName= VcfSpringBeanVariantAnnotator.DEFAULT_MAIN_BEAN_NAME;

	@Override
	protected Logger getLogger() {
		return LOG;
		}
	
	@Override
	protected List<VariantAnnotator> createVariantAnnotators() {
		if(this.springCongigFiles.isEmpty()) throw new IllegalArgumentException("No XML Configuration file was defined");
		final VcfSpringBeanVariantAnnotator ann= new VcfSpringBeanVariantAnnotator(
				springCongigFiles.stream().map(F->F.toUri().toString()).collect(Collectors.toList()),
				mainBeanName
				);
		return Collections.singletonList(ann);
		}
	
	public static void main(final String[] args) {
		new VcfSpringFilter().instanceMainWithExit(args);

	}

}

package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

public class VcfClusteredReadEdgeTest extends TestUtils
	{
	@Test
	public void test01() throws Exception{
		final File output = super.createTmpFile(".vcf");
		Assert.assertEquals(0,new VcfClusteredReadEdge().instanceMain(newCmd().add(
			"-B",SRC_TEST_RESOURCE+"./S1.bam", 
			"-B",SRC_TEST_RESOURCE+"./S2.bam", 
			"-B",SRC_TEST_RESOURCE+"./S3.bam", 
			"-B",SRC_TEST_RESOURCE+"./S4.bam", 
			"-B",SRC_TEST_RESOURCE+"./S5.bam", 
			"-o",output.getPath(),
			SRC_TEST_RESOURCE+"/rotavirus_rf.vcf.gz"
			).make()
			));
		assertIsVcf(output);
		}
	}
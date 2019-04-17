package com.github.lindenb.jvarkit.tools.misc;

import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.vcf.VCFUtilsTest;

@AlsoTest(VCFUtilsTest.class)
public class VcfClusteredReadEdgeTest {
	@Test
	public void test01() throws Exception{
		final TestSupport support = new TestSupport();
		try {
			final Path output = support.createTmpPath(".vcf");
			Assert.assertEquals(new VcfClusteredReadEdge().instanceMain(new String[] {
				"-B",support.resource("S1.bam"), 
				"-B",support.resource("S2.bam"), 
				"-B",support.resource("S3.bam"), 
				"-B",support.resource("S4.bam"), 
				"-B",support.resource("S5.bam"), 
				"-o",output.toString(),
				support.resource("rotavirus_rf.vcf.gz")
				}),0);
			support.assertIsVcf(output);
			} finally {
				support.removeTmpFiles();
			}
		}
	}
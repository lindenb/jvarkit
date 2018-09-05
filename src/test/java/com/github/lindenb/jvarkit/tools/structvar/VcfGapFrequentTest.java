package com.github.lindenb.jvarkit.tools.structvar;

import java.io.File;
import java.io.IOException;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

public class VcfGapFrequentTest extends TestUtils{
	@Test(dataProvider="all-vcf-files")
	public void test01(final String input) 
		throws IOException
		{
		final File out = super.createTmpFile(".bed");
		Assert.assertEquals(new VcfGapFrequent().instanceMain(new String[] {
			"-o",out.getPath(),
			"--fields","AF_NFE",
			"-D",SRC_TEST_RESOURCE+"/gnomad.exomes.r2.0.1.sites.vcf.gz",
			"-D",SRC_TEST_RESOURCE+"/gnomad.genomes.r2.0.1.sites.1.vcf.gz",
			input
			}),0);
		Assert.assertTrue(out.exists());
		}
	}

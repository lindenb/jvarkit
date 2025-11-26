package com.github.lindenb.jvarkit.tools.dict2vcf;

import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;

public class DictToVcfTest {
	
	@Test
	public void basic() throws Exception {
		final TestSupport support = new TestSupport();
		try {
			final Path out = support.createTmpPath(".vcf");
			Assert.assertEquals(new DictToVcf().instanceMain(new String[] {
				"-o",out.toString(),
				support.resource("rotavirus_rf.fa")
				}),0);
			support.assertIsVcf(out);
			}
		finally {
			support.removeTmpFiles();
			}
		}
	
	@Test
	public void withBam() throws Exception {
		final TestSupport support = new TestSupport();
		try {
			final Path out = support.createTmpPath(".vcf");
			Assert.assertEquals(new DictToVcf().instanceMain(new String[] {
				"-o",out.toString(),
				"--samples-file",support.resource("S1.bam"),
				support.resource("rotavirus_rf.fa")
				}),0);
			support.assertIsVcf(out);
			}
		finally {
			support.removeTmpFiles();
			}
		}
}

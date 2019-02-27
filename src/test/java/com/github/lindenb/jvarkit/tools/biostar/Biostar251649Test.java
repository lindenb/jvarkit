package com.github.lindenb.jvarkit.tools.biostar;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;
public class Biostar251649Test {
	
	private final TestSupport support = new TestSupport();

	
	public void test01() throws IOException {
		try {
			final Path out = support.createTmpPath(".vcf");
			Assert.assertEquals(
				new Biostar251649().instanceMain(new String[] {
				"-o",out.toString(),
				"-R",support.resource("rotavirus_rf.fa"),
				support.resource("S1.vcf.gz")
				}),0);
			support.assertIsVcf(out);
			}
		finally
			{
			support.removeTmpFiles();
			}
		}
}

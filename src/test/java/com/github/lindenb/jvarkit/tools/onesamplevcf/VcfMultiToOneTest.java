package com.github.lindenb.jvarkit.tools.onesamplevcf;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.Test;


import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;


public class VcfMultiToOneTest {
	private final TestSupport support = new TestSupport();
	@Test
	public void test01() throws IOException
		{
		try {
			final Path output = support.createTmpPath(".vcf");
			Assert.assertEquals(new VcfMultiToOne().instanceMain(new String[] {
		     		"-o",output.toString(),
		     		support.resource("rotavirus_rf.vcf.gz"),
		     		support.resource("rotavirus_rf.ann.vcf.gz"),
					}),0);
			support.assertIsVcf(output);
			} 
		finally {
			support.removeTmpFiles();
			}
		}
}

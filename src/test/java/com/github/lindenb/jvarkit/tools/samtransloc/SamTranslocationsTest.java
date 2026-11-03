package com.github.lindenb.jvarkit.tools.samtransloc;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;


public class SamTranslocationsTest {

	@Test
	public void simpleTest() throws IOException
		{
		final TestSupport support = new TestSupport();
		try {
			final Path out = support.createTmpPath(".json");
			Assert.assertEquals(new SamTranslocations().instanceMain(new String[] {
				"-o",out.toString(),
				"--chrom-regex","",
				"-R",support.resource("rotavirus_rf.fa"),
				support.resource("S1.bam")//support.resource("HG02260.transloc.chr9.14.bam")
				}),0);
			//support.assertIsVcf(out);
		} finally {
			support.removeTmpFiles();
			}
		}
}

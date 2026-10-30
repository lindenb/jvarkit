package com.github.lindenb.jvarkit.tools.samtransloc;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.samtransloc.SamTranslocations;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;


public class SamTranslocationsTest {

	@Test
	public void simpleTest() throws IOException
		{
		final TestSupport support = new TestSupport();
		try {
			final Path out = support.createTmpPath(".vcf");
			Assert.assertEquals(new SamTranslocations().instanceMain(new String[] {
				"-o",out.toString(),
				support.resource("HG02260.transloc.chr9.14.bam")
				}),0);
			//support.assertIsVcf(out);
		} finally {
			support.removeTmpFiles();
			}
		}
}

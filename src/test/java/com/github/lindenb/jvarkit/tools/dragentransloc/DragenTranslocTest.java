package com.github.lindenb.jvarkit.tools.dragentransloc;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;

public class DragenTranslocTest {
	@Test
	public void simpleTest() throws IOException {
		final TestSupport support = new TestSupport();
		try {
			final Path out = support.createTmpPath(".vcf");
			Assert.assertEquals(new DragenTransloc().instanceMain(new String[] {
				"-o",out.toString(),
				support.resource("dragen.transloc.hg38.vcf")
				}),0);
			support.assertIsVcf(out);
			}
		finally {
			support.removeTmpFiles();
			}
		}
}
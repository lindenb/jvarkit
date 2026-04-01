package com.github.lindenb.jvarkit.tools.vcfwindowsplitter;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.Test;


import com.github.lindenb.jvarkit.tools.tests.TestSupport;



public class VcfSlidingWindowSplitterTest {

	@Test
	public void testZip() throws IOException {
		final TestSupport support = new TestSupport();
		try {
			final Path tmp = support.createTmpPath(".zip");
			final Path manifest = support.createTmpPath(".mf");
			
			Assert.assertEquals(new VcfSlidingWindowSplitter().instanceMain(new String[] {
					"-m",manifest.toString(),
					"-o",tmp.toString(),
					"-w","1000",
					"-s","500",
					"-n","3",
					support.resource("rotavirus_rf.ann.vcf.gz")
					}),0);
			
			support.assertZip(tmp);
			support.assertIsBed(manifest);
			support.assertTsvTableIsConsitent(manifest, null);
			}
		finally {
			support.removeTmpFiles();
		}
	}
}

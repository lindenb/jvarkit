package com.github.lindenb.jvarkit.tools.misc;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;


@AlsoTest(LauncherTest.class)
public class VcfSlidingWindowSplitterTest {
	private final TestSupport support = new TestSupport();

	@Test
	public void test01() throws IOException {
		try {
			final Path tmp = support.createTmpPath(".zip");
			final Path manifest = support.createTmpPath(".mf");
			
			Assert.assertEquals(new VcfGeneSplitter().instanceMain(new String[] {
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

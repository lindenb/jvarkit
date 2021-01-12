package com.github.lindenb.jvarkit.tools.coverage;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;

@AlsoTest(LauncherTest.class)
public class DepthOfCoverageTest {
	private final TestSupport support = new TestSupport();
	
	@Test
	public void test01() throws IOException {
		try {
			final Path out = support.createTmpPath(".tsv");
			Assert.assertEquals(new DepthOfCoverage().instanceMain(new String[] {
				"-o",out.toString(),
				"-R",support.resource("rotavirus_rf.fa"),
				"--auto-mask",
				support.resource("S1.bam"),
				support.resource("S2.bam"),
				}),0);
			support.assertTsvTableIsConsitent(out, null);
		} finally {
			support.removeTmpFiles();
		}
	}
}

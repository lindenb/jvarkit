package com.github.lindenb.jvarkit.tools.bam2graphics;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;

import com.github.lindenb.jvarkit.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;

@AlsoTest(LauncherTest.class)
public class WGSCoveragePlotterTest {
	private final TestSupport support = new TestSupport();
	public void test01() throws IOException {
		try {
			Path out = support.createTmpPath(".svg");

			Assert.assertEquals(new WGSCoveragePlotter().instanceMain(
				new String[]{
				"-R",
				support.resource("rotavirus_rf.fa"),
				"--percentile","median",
				"-C","-1",
				"-o",out.toString(),
				"--partition","sample",
				"--samples","S1,Sx",
				"--points",
				support.resource("S1.bam")
				}),
				0);
			support.assertIsXml(out);
			}
		finally
			{
			support.removeTmpFiles();
			}
		}

}

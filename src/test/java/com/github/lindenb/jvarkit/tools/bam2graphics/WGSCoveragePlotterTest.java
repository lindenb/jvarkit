package com.github.lindenb.jvarkit.tools.bam2graphics;

import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.testng.Assert;

import com.github.lindenb.jvarkit.tools.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;

@AlsoTest(Launcher.class)
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

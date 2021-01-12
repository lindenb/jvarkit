package com.github.lindenb.jvarkit.tools.bam2graphics;

import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.testng.Assert;

import com.github.lindenb.jvarkit.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;

@AlsoTest(LauncherTest.class)
public class CoveragePlotterTest {
	private final TestSupport support = new TestSupport();
	public void testZip() throws IOException {
		try {
			final Path imgOut = support.createTmpPath(".zip");
			final List<String> args= new ArrayList<>();
			args.add("-R");
			args.add(support.resource("rotavirus_rf.fa"));
			args.add("-o");
			args.add(imgOut.toString());
			Arrays.asList("1","2","3","4","5").stream().
				map(S->support.resource("/S"+S+".bam")).
				forEach(S->{
					args.add("-B");
					args.add(S);
				});
			args.add("RF01:1-1000");

			Assert.assertEquals(new CoveragePlotter().instanceMain(
				args.toArray(new String[args.size()])),
				0);
			support.assertZip(imgOut);
			}
		finally
			{
			support.removeTmpFiles();
			}
		}

}

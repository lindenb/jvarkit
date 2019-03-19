package com.github.lindenb.jvarkit.tools.structvar;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;


@AlsoTest(LauncherTest.class)
public class SamScanSplitReadsTest {
	@Test
	public void test01() throws java.io.IOException {
		final TestSupport support = new TestSupport();

		try {
		final Path out= support.createTmpPath(".vcf");
		Assert.assertEquals(new SamScanSplitReads().instanceMain(new String[] {
				"-o",out.toString(),
				support.resource("S1.bam"),
				support.resource("S2.bam"),
				support.resource("S3.bam")
				}),0);
		support.assertIsVcf(out);
		} finally
		{
			support.removeTmpFiles();
		}
	}

}

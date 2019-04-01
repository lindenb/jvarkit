package com.github.lindenb.jvarkit.tools.structvar;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;

@AlsoTest(LauncherTest.class)
public class SamTranslocationsTest {
	final TestSupport support = new TestSupport();

	@Test
	public void test01() 
		throws IOException
		{
		try {
			final Path out = support.createTmpPath(".vcf");
			Assert.assertEquals(new SamTranslocations().instanceMain(new String[] {
				"-o",out.toString(),
				support.resource("HG02260.transloc.chr9.14.bam")
				}),0);
			support.assertIsVcf(out);
		} finally {
			support.removeTmpFiles();
		}
		}
}

package com.github.lindenb.jvarkit.tools.cmpbams;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;

@AlsoTest(LauncherTest.class)
public class BamMatrixTest {
	private final TestSupport support = new TestSupport();
	@Test
	public void test01() throws IOException
		{
		try {
		final Path out = support.createTmpPath(".png");
		Assert.assertEquals(new BamMatrix().instanceMain(new String[] {
			"-o",out.toString(),
			"-r","RF01",
			support.resource("S1.bam")
			}),0);
		Assert.assertTrue(support.isImage(out));
		} finally {
			support.removeTmpFiles();
			}
		}

}

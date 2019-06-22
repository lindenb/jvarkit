package com.github.lindenb.jvarkit.tools.cmpbams;

import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;

public class BamMatrixTest {
	private final TestSupport support = new TestSupport();
	@Test(dataProvider="src1")
	public void test01() throws IOException
		{
		try {
		final Path out = support.createTmpPath(".png");
		Assert.assertEquals(new BamMatrix().instanceMain(new String[] {
			"-o",out.toString(),
			"-r","RF01",
			support.resource("S1.bam")
			}),0);
		} finally {
			support.removeTmpFiles();
			}
		}

}

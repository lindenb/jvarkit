package com.github.lindenb.jvarkit.tools.bammatrix;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.Test;


import com.github.lindenb.jvarkit.tools.tests.TestSupport;


public class BamMatrixTest {
	@Test
	public void bamMatrixSimple() throws IOException
		{
		final TestSupport support = new TestSupport();

		try {
		final Path out = support.createTmpPath(".png");
		Assert.assertEquals(new BamMatrix().instanceMain(new String[] {
			"-o",out.toString(),
			"-r","RF01:1-1000",
			support.resource("S1.bam")
			}),0);
		Assert.assertTrue(support.isImage(out));
		} finally {
			support.removeTmpFiles();
			}
		}

}

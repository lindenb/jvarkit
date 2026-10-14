package com.github.lindenb.jvarkit.tools.bed2svg;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;

public class BedToSVGTest {
	@Test
	public void testWithDict()
		throws IOException
		{
		final TestSupport support = new TestSupport();
		try {
			final Path out = support.createTmpPath(".svg");
			final BedToSVG cmd =new BedToSVG();
			Assert.assertEquals(0,cmd.instanceMain(new String[] {
				"-R",support.resource("rotavirus_rf.dict"),
				"-o",out.toString(),
				support.resource("rotavirus_rf.vcf.gz"),
				}));
			support.assertIsXml(out);
			}
		finally {
			support.removeTmpFiles();
			}
		}
	}
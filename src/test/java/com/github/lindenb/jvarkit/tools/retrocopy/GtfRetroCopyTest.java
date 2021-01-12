package com.github.lindenb.jvarkit.tools.retrocopy;


import com.github.lindenb.jvarkit.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.bio.structure.GtfReaderTest;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.Test;

@AlsoTest({LauncherTest.class,GtfReaderTest.class})
public class GtfRetroCopyTest   {
	private final TestSupport support = new TestSupport();
	
	@Test
	public void test01() 
			throws IOException
			{
			try {
				final Path out = support.createTmpPath(".vcf");
				Assert.assertEquals(new GtfRetroCopy().instanceMain(new String[] {
						"--gtf",support.resource("Homo_sapiens.GRCh37.87.gtf.gz"),
						"-o",out.toString(),
						support.resource("test_vcf01.vcf")
						}),0);
				
				support.assertIsVcf(out);
				}
			finally {
				support.removeTmpFiles();
				}
			}
}

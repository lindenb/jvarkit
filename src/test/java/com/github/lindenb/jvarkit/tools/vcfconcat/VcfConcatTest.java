package com.github.lindenb.jvarkit.tools.vcfconcat;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;

public class VcfConcatTest {


	@Test
	public void testMerge() throws IOException
		{
		final TestSupport support=new TestSupport();
		try {
			final Path output = support.createTmpPath(".vcf");
			Assert.assertEquals(new VcfConcat().instanceMain(new String[] {
	        		"-o",output.toString(),
	        		"--merge",
	        		 support.resource("S1.vcf.gz"),
	        		 support.resource("S2.vcf.gz"),
	        		 support.resource("S3.vcf.gz"),
	        		 support.resource("S4.vcf.gz"),
	        		 support.resource("S5.vcf.gz")
				}),0);
			support.assertIsVcf(output);
			}
		finally {
			support.removeTmpFiles();
			}
		}
	}

package com.github.lindenb.jvarkit.tools.misc;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.tools.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.vcf.VCFUtilsTest;

@AlsoTest(VCFUtilsTest.class)
public class NoZeroVariationVCFTest {
	
	@Test
	public void test01() 
		throws IOException
		{
		final TestSupport support = new TestSupport();
		try {
			final Path inputFile = support.createTmpPath(".vcf");
			IOUtils.cat("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n",
					inputFile, false);
			
			final Path out = support.createTmpPath(".vcf");
			Assert.assertEquals(0,new NoZeroVariationVCF().instanceMain(new String[] {
				"-o",out.toString(),
				"-R",support.resource("rotavirus_rf.fa"),
				inputFile.toString()
				}));
			support.assertIsVcf(out);
		} finally {
			support.removeTmpFiles();
		}
		}
}

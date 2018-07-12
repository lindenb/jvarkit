package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

import org.testng.Assert;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

public class NoZeroVariationVCFTest extends TestUtils {
	
	public void test01() 
		throws IOException
		{
		final File inputFile = super.createTmpFile(".vcf");
		PrintWriter pw= new PrintWriter(inputFile);
		pw.println("##fileformat=VCFv4.2");
		pw.println("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1");
		pw.flush();
		pw.close();
		
		final File out = super.createTmpFile(".vcf");
		Assert.assertEquals(0,new NoZeroVariationVCF().instanceMain(new String[] {
			"-o",out.getPath(),
			"-R",SRC_TEST_RESOURCE+"/rotavirus_rf.fa",
			inputFile.getPath()
			}));
		assertIsVcf(out);
		}
}

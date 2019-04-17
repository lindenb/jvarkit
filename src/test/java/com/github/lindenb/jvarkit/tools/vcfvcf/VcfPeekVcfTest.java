package com.github.lindenb.jvarkit.tools.vcfvcf;

import java.io.IOException;
import java.nio.file.Path;
import java.util.Arrays;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;

@AlsoTest(LauncherTest.class)
public class VcfPeekVcfTest {
	private final TestSupport support =new TestSupport();

@DataProvider(name="src01")
public Object[][] src01(){
	return new Object[][] {
		{support.resource("S1.vcf.gz"),support.resource("S2.vcf.gz")},
		{support.resource("S2.vcf.gz"),support.resource("S3.vcf.gz")},
		{support.resource("S3.vcf.gz"),support.resource("S4.vcf.gz")},
		{support.resource("S4.vcf.gz"),support.resource("S5.vcf.gz")},
		{support.resource("S5.vcf.gz"),support.resource("S1.vcf.gz")}
	};
}
	
@Test(dataProvider="src01")
public void test01(final String vcfIn,final String vcfdb)
	throws IOException
	{
	try {
		final Path out = support.createTmpPath(".vcf"); 
		final String args[]=new String[] {
				"-o",out.toString(),
				"-f",vcfdb,
				"-t","AN,AC,DP",
				"-p","TITITOTO",
				"--alt","at_least_one",
				vcfIn
				};
		int ret = new VcfPeekVcf().instanceMain(args);
		if(ret!=0)
			{
			Assert.fail(Arrays.asList(args).toString());
			}
		support.assertIsVcf(out);
		}
	finally {
		support.removeTmpFiles();
		}
	}
	
}

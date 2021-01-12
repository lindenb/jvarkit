package com.github.lindenb.jvarkit.tools.vcffilterjs;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;

@AlsoTest(LauncherTest.class)
public class VcfFilterJdkTest {
private final TestSupport support =new TestSupport();
	
	@DataProvider(name = "src1")
	public Object[][] createData1() {
		return new Object[][] {
			{support.resource("ExAC.r1.sites.vep.vcf.gz"),"return variant.getStart()%10==0;"},
			{support.resource("rotavirus_rf.ann.vcf.gz"),"return variant.getContig().equals(\"RF01\");"}
		};
		}

	
	@Test(dataProvider="src1")
	public void test01(final String inputFile,String expr) 
		throws IOException
		{
		try {
			final Path output = support.createTmpPath(".vcf");
	        Assert.assertEquals(new VcfFilterJdk().instanceMain(new String[] {
	        		"-o",output.toString(),
	        		"-e",expr,
	        		inputFile}),
	        		0);
	        support.assertIsVcf(output);
			} 
		finally
			{	
			support.removeTmpFiles();
			}
		}
}
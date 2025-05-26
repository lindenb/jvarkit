package com.github.lindenb.jvarkit.tools.vcffilterjs;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;


import com.github.lindenb.jvarkit.tools.tests.TestSupport;


public class VcfFilterJdkTest {

	
	@DataProvider(name = "src1")
	public Object[][] createData1() {
		final TestSupport support =new TestSupport();
		
		return support.allVcfOrBcf().
				filter(F->F.endsWith(".vcf") || F.endsWith(".vcf.gz")).
				flatMap(F->Arrays.stream(new Object[][]{
					{F,"return variant.getStart()%10==0;"},
					{F,"return variant.getContig().equals(\"RF01\");"}
					})).toArray(N->new Object[N][]);
		}

	
	@Test(dataProvider="src1")
	public void testExpression(final String inputFile,String expr) 
		throws IOException
		{
		final TestSupport support =new TestSupport();
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
	
	@Test(dataProvider="src1")
	public void testFile(final String inputFile,String expr) 
		throws IOException
		{
		final TestSupport support =new TestSupport();
		try {
			Path code  = support.createTmpPath(".code");
			Files.writeString(code, expr);
			final Path output = support.createTmpPath(".vcf");
	        Assert.assertEquals(new VcfFilterJdk().instanceMain(new String[] {
	        		"-o",output.toString(),
	        		"-f",code.toAbsolutePath().toString(),
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
package com.github.lindenb.jvarkit.tools.vcffilterso;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;

@AlsoTest(LauncherTest.class)
public class VcfFilterSequenceOntologyTest {
	private final TestSupport support =new TestSupport();
	
	@DataProvider(name = "src1")
	public Object[][] createData1() {
		return new Object[][] {
			{support.resource("ExAC.r1.sites.vep.vcf.gz")},
			{support.resource("rotavirus_rf.ann.vcf.gz")}
		};
		}

	
	@Test(dataProvider="src1")
	public void test01(final String inputFile) 
		throws IOException
		{
		try {

			final Path output = support.createTmpPath(".vcf");
	        Assert.assertEquals(new VcfFilterSequenceOntology().instanceMain(new String[] {
	        		"-o",output.toString(),
	        		"-A","SO:0001818",
	        		inputFile}),
	        		0);
	        support.assertIsVcf(output);
	        Assert.assertTrue(support.wc(output)<=support.wc(Paths.get(inputFile)));
			} 
		finally
			{	
			support.removeTmpFiles();
			}
		}
	
	@Test(dataProvider="src1")
	public void test02(final String inputFile) 
		throws IOException
		{
		try {

			final Path output = support.createTmpPath(".vcf");
	        Assert.assertEquals(new VcfFilterSequenceOntology().instanceMain(new String[] {
	        		"-o",output.toString(),
	        		"-A","SO:0001818",
	        		"--filterin","XXX",
	        		inputFile}),
	        		0);
	        support.assertIsVcf(output);
	        Assert.assertEquals(support.wc(output),support.wc(Paths.get(inputFile)));
			} 
		finally
			{	
			support.removeTmpFiles();
			}
		}
	
	@Test(dataProvider="src1")
	public void test03(final String inputFile) 
		throws IOException
		{
		try {

			final Path output = support.createTmpPath(".vcf");
	        Assert.assertEquals(new VcfFilterSequenceOntology().instanceMain(new String[] {
	        		"-o",output.toString(),
	        		"-A","SO:0001818",
	        		"--filterout","XXX",
	        		inputFile}),
	        		0);
	        support.assertIsVcf(output);
	        Assert.assertEquals(support.wc(output),support.wc(Paths.get(inputFile)));
			} 
		finally
			{	
			support.removeTmpFiles();
			}
		}

}

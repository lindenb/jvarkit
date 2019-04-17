package com.github.lindenb.jvarkit.tools.biostar;

import java.nio.file.Path;
import java.nio.file.Paths;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;

public class Biostar322664Test {
	
	private final TestSupport support = new TestSupport();

	
	@DataProvider(name = "src01")
	public Object[][] getData01() {
		return new Object[][]{
			{support.resource("S1.bam"),support.resource("S1.vcf.gz")},
			{support.resource("S2.bam"),support.resource("S2.vcf.gz")},
			{support.resource("S3.bam"),support.resource("S3.vcf.gz")}
			};
		}
	
	@Test(dataProvider="src01")
	public void testWithQueryNameSortedBam(final String bampath,final String vcfpath) throws Exception
		{
			try {
			final Path out = support.createTmpPath(".bam");
			final Path sortedBam1= support.sortBamOnQueryName(Paths.get(bampath),null);
			
			Assert.assertEquals(new Biostar322664().instanceMain(new String[] {
				"-o",out.toString(),
				"-V",vcfpath,
				sortedBam1.toString()
				}),0);
			support.assertIsValidBam(out);
			} 
		finally
			{
			support.removeTmpFiles();
			}
		}
	
	@Test(dataProvider="src01")
	public void testWithBamIndex(final String bampath,final String vcfpath) throws Exception
		{
		try {
			final Path out = support.createTmpPath(".bam");
			
			Assert.assertEquals(new Biostar322664().instanceMain(new String[] {
				"-o",out.toString(),"-index","-nm",
				"-V",vcfpath,
				bampath
				}),0);
			support.assertIsValidBam(out);
			} 
		finally
			{
			support.removeTmpFiles();
			}
		}
	}

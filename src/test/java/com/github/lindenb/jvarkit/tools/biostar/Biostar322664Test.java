package com.github.lindenb.jvarkit.tools.biostar;

import java.io.File;
import java.nio.file.Paths;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

public class Biostar322664Test extends TestUtils{
	@DataProvider(name = "src01")
	public Object[][] getData01() {
		return new Object[][]{
			{SRC_TEST_RESOURCE+"/S1.bam",SRC_TEST_RESOURCE+"/S1.vcf.gz"},
			{SRC_TEST_RESOURCE+"/S2.bam",SRC_TEST_RESOURCE+"/S2.vcf.gz"},
			{SRC_TEST_RESOURCE+"/S3.bam",SRC_TEST_RESOURCE+"/S3.vcf.gz"}
			};
		}
	
	@Test(dataProvider="src01")
	public void testWithQueryNameSortedBam(final String bampath,final String vcfpath) throws Exception
		{
		final File out = super.createTmpFile(".bam");
		final File sortedBam1= sortBamOnQueryName(Paths.get(bampath),null);
		
		Assert.assertEquals(new Biostar322664().instanceMain(new String[] {
			"-o",out.getPath(),
			"-V",vcfpath,
			sortedBam1.getPath()
			}),0);
		super.assertIsValidBam(out);
		}
	
	@Test(dataProvider="src01")
	public void testWithBamIndex(final String bampath,final String vcfpath) throws Exception
		{
		final File out = super.createTmpFile(".bam");
		
		Assert.assertEquals(new Biostar322664().instanceMain(new String[] {
			"-o",out.getPath(),"-index","-nm",
			"-V",vcfpath,
			bampath
			}),0);
		super.assertIsValidBam(out);
		}
	}

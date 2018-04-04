package com.github.lindenb.jvarkit.tools.vcfvcf;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

public class VcfPeekVcfTest extends TestUtils{

@DataProvider(name="src01")
public Object[][] src01(){
	return new Object[][] {
		{SRC_TEST_RESOURCE+"/S1.vcf.gz",SRC_TEST_RESOURCE+"/S2.vcf.gz"},
		{SRC_TEST_RESOURCE+"/S2.vcf.gz",SRC_TEST_RESOURCE+"/S3.vcf.gz"},
		{SRC_TEST_RESOURCE+"/S3.vcf.gz",SRC_TEST_RESOURCE+"/S4.vcf.gz"},
		{SRC_TEST_RESOURCE+"/S4.vcf.gz",SRC_TEST_RESOURCE+"/S5.vcf.gz"},
		{SRC_TEST_RESOURCE+"/S5.vcf.gz",SRC_TEST_RESOURCE+"/S1.vcf.gz"}
	};
}
	
@Test(dataProvider="src01")
public void test01(final String vcfIn,final String vcfdb)
	throws IOException
	{
	final File out = super.createTmpFile(".vcf"); 
	final String args[]=newCmd().
			add("-o",out).
			add("-f",vcfdb).
			add("-t","AN,AC,DP").
			add("-p","TITITOTO").
			add("--alt","at_least_one").
			add(vcfIn).
			make();
	int ret = new VcfPeekVcf().instanceMain(args);
	if(ret!=0)
		{
		Assert.fail(Arrays.asList(args).toString());
		}
	assertIsVcf(out);
	}
	
}

package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.IOException;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

public class PadEmptyFastqTest extends TestUtils {
	
	@DataProvider(name="src01")
	public Object[][] getData() {
		return new Object[][] {
			{SRC_TEST_RESOURCE+"/S1.R1.fq.gz"},
			{SRC_TEST_RESOURCE+"/S1.R2.fq.gz"},
			{SRC_TEST_RESOURCE+"/S2.R1.fq.gz"},
			{SRC_TEST_RESOURCE+"/S2.R2.fq.gz"},
			{SRC_TEST_RESOURCE+"/S3.R1.fq.gz"},
			{SRC_TEST_RESOURCE+"/S3.R2.fq.gz"},
			{SRC_TEST_RESOURCE+"/S4.R1.fq.gz"},
			{SRC_TEST_RESOURCE+"/S4.R2.fq.gz"},
			{SRC_TEST_RESOURCE+"/S5.R1.fq.gz"},
			{SRC_TEST_RESOURCE+"/S5.R2.fq.gz"}
		};
	}
	
	@Test(dataProvider="src01")
	public void test01(final String fq) 
			throws IOException
			{
			final File out = super.createTmpFile(".fq");
		
			Assert.assertEquals(0,new PadEmptyFastq().instanceMain(new String[] {
				"-o",out.getPath(),
				fq
				}));
			assertIsFastq(out);
			}
}

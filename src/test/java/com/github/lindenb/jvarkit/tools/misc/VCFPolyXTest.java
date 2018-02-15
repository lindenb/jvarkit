package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.IOException;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

public class VCFPolyXTest extends TestUtils{
		@DataProvider(name = "src1")
		public Object[][] createData1() {
			return new Object[][]{
				{"./src/test/resources/toy.vcf.gz","./src/test/resources/toy.fa",1}
				};
		}
	
		
	@Test(dataProvider="src1")
	public void test1(final String vcf,final String ref,int n) throws IOException {
		final File out = createTmpFile(".vcf");
		Assert.assertEquals(
			new VCFPolyX().instanceMain(newCmd().add(
					"-o",out,
	        		"-R",ref,
	        		"-n","1",
	        		vcf
			).make()),0);
		assertIsVcf(out);
		}
}

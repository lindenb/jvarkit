package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.IOException;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

public class VcfAlleleBalanceTest extends TestUtils {
	@Test(dataProvider="all-vcf-files")
	public void test01(final String inputFile) 
		throws IOException
		{
		final File out = super.createTmpFile(".vcf");
		Assert.assertEquals(0,new VcfAlleleBalance().instanceMain(new String[] {
			"-o",out.getPath(),
			inputFile
			}));
		assertIsVcf(out);
		}
	@Test(dataProvider="all-vcf-files")
	public void testWithPed(final String inputFile) 
		throws IOException
		{
		final File out = super.createTmpFile(".vcf");
		final File ped = super.createRandomPedigreeFromFile(inputFile);
		if(ped==null) return;
		Assert.assertEquals(0,new VcfAlleleBalance().instanceMain(new String[] {
			"-o",out.getPath(),
			"-p",ped.getPath(),
			inputFile
			}));
		assertIsVcf(out);
		}
	}

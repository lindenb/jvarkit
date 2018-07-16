package com.github.lindenb.jvarkit.tools.vcffixindels;

import java.io.File;
import java.io.IOException;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

public class VCFFixIndelsTest extends TestUtils {
	@Test(dataProvider="all-vcf-files")
	public void test01(final String inputFile) 
		throws IOException
		{
		final File output = super.createTmpFile(".vcf");
		Assert.assertEquals(new VCFFixIndels().instanceMain(
        		newCmd().add(
        		"-o",output,
        		inputFile).make()
        	),0);
        assertIsVcf(output);
		}
	}

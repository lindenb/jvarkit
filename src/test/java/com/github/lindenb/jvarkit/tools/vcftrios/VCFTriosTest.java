package com.github.lindenb.jvarkit.tools.vcftrios;

import java.io.File;
import java.io.IOException;

import org.testng.Assert;
import org.testng.Reporter;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

public class VCFTriosTest extends TestUtils {
	@Test(dataProvider="all-vcf-files")
	public void test01(final String inputFile) 
		throws IOException
		{
		final File ped = super.createRandomPedigreeFromFile(inputFile);
		if(ped==null) {
			 Reporter.getCurrentTestResult().setAttribute("warn", "No Pedigree for "+inputFile);
			return;
		}
		final File output = super.createTmpFile(".vcf");
		
		Assert.assertEquals(new VCFTrios().instanceMain(
        		newCmd().add(
        		"-o",output,
        		"--pedigree",ped,
        		inputFile).make()
        	),0);
        assertIsVcf(output);
		}
	}

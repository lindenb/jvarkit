package com.github.lindenb.jvarkit.tools.burden;

import java.io.File;
import java.io.IOException;

import org.testng.Assert;
import org.testng.Reporter;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

public class VcfFilterNotInPedigreeTest extends TestUtils {

	@Test(dataProvider="all-vcf-files")
	public void test01(final String inputFile) 
		throws IOException
		{
		final VCFFileReader r0 = new VCFFileReader(new File(inputFile),false);
		final VCFHeader vcfheader = r0.getFileHeader();
		if(vcfheader.getNGenotypeSamples()<1) {
			r0.close();
			return;
		}
		r0.close();
		
		final File ped = super.createRandomPedigreeFromFile(inputFile);
		if(ped==null) {
			 Reporter.getCurrentTestResult().setAttribute("warn", "No Pedigree for "+inputFile);
			return;
		}
		final File output = super.createTmpFile(".vcf");
		
		Assert.assertEquals(new VcfFilterNotInPedigree().instanceMain(
        		newCmd().add(
        		"-o",output,
        		"--pedigree",ped,
        		inputFile).make()
        	),0);
        assertIsVcf(output);
		}


}

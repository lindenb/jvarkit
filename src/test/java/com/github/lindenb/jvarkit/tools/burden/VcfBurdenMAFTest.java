package com.github.lindenb.jvarkit.tools.burden;

import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;

import org.testng.Assert;
import org.testng.Reporter;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;

import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

@AlsoTest(LauncherTest.class)
public class VcfBurdenMAFTest {
	
	private final TestSupport support = new TestSupport();

	
	@Test
	public void test01() 
		throws IOException
		{
		try {
			final String inputFile = support.resource("rotavirus_rf.vcf.gz");
			
			try(final VCFFileReader r0 = new VCFFileReader(Paths.get(inputFile),false)) {
					final VCFHeader vcfheader = r0.getFileHeader();
					if(vcfheader.getNGenotypeSamples()<1) {
						return;
					}
				}
			
			final Path ped = support.createRandomPedigreeFromFile(inputFile);
			if(ped==null) {
				 Reporter.getCurrentTestResult().setAttribute("warn", "No Pedigree for "+inputFile);
				return;
			}
			final Path output = support.createTmpPath(".vcf");
			
			Assert.assertEquals(new VcfBurdenMAF().instanceMain(new String[] {
	        		"-o",output.toString(),
	        		"--max-af","5/100",
	        		"--pedigree",ped.toString(),
	        		inputFile.toString()}),0);
			support.assertIsVcf(output);
			} 
		finally {
			support.removeTmpFiles();
			}
		}
}

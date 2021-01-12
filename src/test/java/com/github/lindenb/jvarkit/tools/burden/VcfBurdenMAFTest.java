package com.github.lindenb.jvarkit.tools.burden;

import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;

import org.testng.Assert;
import org.testng.Reporter;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;
import com.github.lindenb.jvarkit.variant.vcf.VCFReaderFactory;

import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFReader;

@AlsoTest(LauncherTest.class)
public class VcfBurdenMAFTest {
	
	private final TestSupport support = new TestSupport();

	
	@Test
	public void test01() 
		throws IOException
		{
		try {
			final String inputFile = support.resource("rotavirus_rf.vcf.gz");
			
			try(final VCFReader r0 = VCFReaderFactory.makeDefault().open(Paths.get(inputFile),false)) {
					final VCFHeader vcfheader = r0.getHeader();
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

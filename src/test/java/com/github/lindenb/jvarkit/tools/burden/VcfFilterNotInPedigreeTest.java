package com.github.lindenb.jvarkit.tools.burden;

import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;

import org.testng.Assert;
import org.testng.Reporter;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;

import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

@AlsoTest(LauncherTest.class)
public class VcfFilterNotInPedigreeTest  {
	private final TestSupport support = new TestSupport();

	@DataProvider(name = "src1")
	public Object[][] createData1() {
		return (Object[][])support.
				allVcfOrBcf().
				map(F->new Object[] {F}).
				toArray()
				;
		} 
	
	@Test(dataProvider="src1")
	public void test01(final String inputFile) 
		throws IOException
		{
		try {
			final VCFFileReader r0 = new VCFFileReader(Paths.get(inputFile),false);
			final VCFHeader vcfheader = r0.getFileHeader();
			if(vcfheader.getNGenotypeSamples()<1) {
				r0.close();
				return;
			}
			r0.close();
			
			final Path ped = support.createRandomPedigreeFromFile(inputFile);
			if(ped==null) {
				 Reporter.getCurrentTestResult().setAttribute("warn", "No Pedigree for "+inputFile);
				return;
			}
			final Path output = support.createTmpPath(".vcf");
			
			Assert.assertEquals(new VcfFilterNotInPedigree().instanceMain(new String[] {
	        		"-o",output.toString(),
	        		"--pedigree",ped.toString(),
	        		inputFile.toString()}),0);
			support.assertIsVcf(output);
			} 
		finally {
			support.removeTmpFiles();
			}
		}


}

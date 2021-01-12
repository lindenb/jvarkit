package com.github.lindenb.jvarkit.tools.burden;

import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;

import org.testng.Assert;
import org.testng.Reporter;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;
import com.github.lindenb.jvarkit.variant.vcf.VCFReaderFactory;

import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFReader;

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
			final VCFReader r0 =  VCFReaderFactory.makeDefault().open(Paths.get(inputFile),false);
			final VCFHeader vcfheader = r0.getHeader();
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

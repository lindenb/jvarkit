package com.github.lindenb.jvarkit.tools.burden;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;


import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.tools.vcf2r.VcfToRScript;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;


public class VcfBurdenRscriptVTest  {
	
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
			final Path ped = support.createRandomPedigreeFromFile(inputFile);
			if(ped==null) {
				return;
				}
			final Path output = support.createTmpPath(".vcf");
			
			Assert.assertEquals(new VcfToRScript().instanceMain(new String[] {
	        		"-o",output.toString(),
	        		"--pedigree",ped.toString(),
					}),0);
			support.assertIsNotEmpty(output);
			} 
		finally {
			support.removeTmpFiles();
			}
		}

}

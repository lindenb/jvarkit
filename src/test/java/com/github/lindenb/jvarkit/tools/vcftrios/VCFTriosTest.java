package com.github.lindenb.jvarkit.tools.vcftrios;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.Reporter;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;

@AlsoTest({DeNovoDetector.class,LauncherTest.class})
public class VCFTriosTest {
	
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
	public void testDeNovo(final String inputFile) 
		throws IOException
		{
		try {
		final Path ped = support.createRandomPedigreeFromFile(inputFile);
			if(ped==null) {
				 Reporter.getCurrentTestResult().setAttribute("warn", "No Pedigree for "+inputFile);
				return;
			}
			final Path output = support.createTmpPath(".vcf");
			
			Assert.assertEquals(new VCFTrios().instanceMain(new String[] {
	        		"-o",output.toString(),
	        		"--pedigree",ped.toString(),
	        		inputFile}),0);
			support.assertIsVcf(output);
			} 
		finally {
			support.removeTmpFiles();
			}
		}
	}

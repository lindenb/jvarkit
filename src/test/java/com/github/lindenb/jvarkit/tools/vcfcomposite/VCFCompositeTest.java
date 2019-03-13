package com.github.lindenb.jvarkit.tools.vcfcomposite;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;

@AlsoTest(LauncherTest.class)
public class VCFCompositeTest {
	private final TestSupport support =new TestSupport();

	@DataProvider(name = "src1")
	public Object[][] createData1() {
		return support.toArrayArray(support.
				allVcfOrBcf().
				map(F->new Object[] {F})
				)
				;
		}
	
	@Test(dataProvider="src1")
	public void test01(final String inputFile) 
		throws IOException
		{
		try {
			final Path ped = support.createRandomPedigreeFromFile(inputFile);
			if(ped==null) return;
			if(Files.lines(ped).noneMatch(L->L.endsWith("1")))
				{
				//no affected in pedigree ?
				return;
				}
			if(Files.lines(ped).noneMatch(L->L.endsWith("0")))
				{
				//no unaffected in pedigree ?
				return;
				}
			final Path output = support.createTmpPath(".vcf");
	        Assert.assertEquals(new VCFComposite().instanceMain(new String[] {
	        		"-o",output.toString(),
	        		"--pedigree",ped.toString(),
	        		inputFile}),0);
	        support.assertIsVcf(output);
			} 
		finally
			{	
			support.removeTmpFiles();
			}
		}
}

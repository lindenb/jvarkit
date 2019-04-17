package com.github.lindenb.jvarkit.tools.bioalcidae;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;

public class BioAlcidaeTest
	{
	private final TestSupport support = new TestSupport();

	@DataProvider(name = "src1")
	public Object[][] createData1() {
		return (Object[][])support.
				allVcfOrBcf().
				map(F->new Object[] {F}).
				toArray()
				;
		}
	
	@Test(dataProvider = "src1")
	public void testVCF(final String vcfpath) 
		throws IOException
		{
		try {
			final Path out = support.createTmpPath(".txt");
			
			Assert.assertEquals(new BioAlcidae().instanceMain(new String[] {
					"-o",
					out.toString(),
					"-e",
					"while(iter.hasNext()) out.println(iter.next().getContig());",
					vcfpath
					}),0);
			support.assertIsNotEmpty(out);
			} 
		finally
			{
			support.removeTmpFiles();
			}
		}
	
}

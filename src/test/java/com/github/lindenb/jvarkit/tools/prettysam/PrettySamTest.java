package com.github.lindenb.jvarkit.tools.prettysam;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;


public class PrettySamTest {
	@DataProvider(name = "src1")
	public Object[][] createData1() {
		final TestSupport support = new TestSupport();
		return support.toArrayArray(support.
				allSamOrBams().
				map(F->new Object[] {F})
				)
				;
		}
	
	@Test(dataProvider="src1")
	public void prettySam1(final String inBam,String inFasta) 
		throws IOException
		{
		final TestSupport support = new TestSupport();
		try {
			final Path out = support.createTmpPath(".txt");
			final PrettySam cmd =new PrettySam();
			Assert.assertEquals(cmd.instanceMain(new String[] {
				"-R",inFasta,
				"-o",out.toString(),
				inBam
				}),0);
			support.assertIsNotEmpty(out);
			}
		finally {
			support.removeTmpFiles();
			}
		}
}

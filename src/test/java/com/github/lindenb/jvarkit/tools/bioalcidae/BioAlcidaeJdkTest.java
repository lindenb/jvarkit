package com.github.lindenb.jvarkit.tools.bioalcidae;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;

public class BioAlcidaeJdkTest
	{
	private final TestSupport support = new TestSupport();

	@DataProvider(name = "src1")
	public Object[][] createData1() {
		return (Object[][])support.
				allVcfOrBcf().
				map(F->new Object[] {F,"stream().forEach(V->out.println(V.getContig()));"}).
				toArray()
				;
		}
	
	private void simpleTest(final String inputFile,final String script) throws IOException
		{
		try {
			final Path out = support.createTmpPath(".txt");
			Assert.assertEquals(0,new BioAlcidaeJdk().instanceMain(new String[] {
				"-e",script,
				"-o",out.toString(),
				inputFile
				}),0);
			support.assertIsNotEmpty(out);
			}
		finally
			{
			support.removeTmpFiles();
			}
		}

	@Test(dataProvider="src1")
	public void testVcf(final String inputFile,final String script) 
		throws IOException
		{
		simpleTest(inputFile,script);
		}

	
	@DataProvider(name = "src2")
	public Object[][] createData2() {
		return (Object[][])support.
				allSamOrBams().
				map(F->new Object[] {F,"stream().forEach(V->out.println(V.getContig()));"}).
				toArray()
				;
		}

	@Test(dataProvider="src2")
	public void testBam(final String inputFile,final String script) 
		throws IOException
		{
		simpleTest(inputFile,script);
		}

	@DataProvider(name = "src3")
	public Object[][] createData3() {
		return (Object[][])support.
				allFasta().
				map(F->new Object[] {F,"stream().forEach(V->out.println(V.length()));"}).
				toArray()
				;
		}

	@Test(dataProvider="src3")
	public void testFasta(final String inputFile,final String script) 
		throws IOException
		{
		simpleTest(inputFile,script);
		}
	
}

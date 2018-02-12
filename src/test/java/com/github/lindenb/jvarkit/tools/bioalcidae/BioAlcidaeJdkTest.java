package com.github.lindenb.jvarkit.tools.bioalcidae;

import java.io.File;
import java.io.IOException;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

public class BioAlcidaeJdkTest extends TestUtils
	{
	@DataProvider(name = "src1")
	public Object[][] createData1() {
		return new ParamCombiner().
			initList(collectAllVcfs()).
			product(new Object[] {
					"stream().forEach(V->out.println(V.getContig()));"
			}).
			build();
	}
	
	private void simpleTest(final String inputFile,final String script) throws IOException
		{
		final File out = File.createTempFile(".tmp.", ".txt");
		final BioAlcidaeJdk cmd =new BioAlcidaeJdk();
		Assert.assertEquals(0,cmd.instanceMain(new String[] {
			"-e",script,
			"-o",out.getPath(),
			inputFile
			}));
		Assert.assertTrue(out.delete());
		}

	@Test(dataProvider="src1")
	public void test01(final String inputFile,final String script) 
		throws IOException
		{
		simpleTest(inputFile,script);
		}

	
	@DataProvider(name = "src2")
	public Object[][] createData2() {
		return new ParamCombiner().
			initList(collectAllSamOrBam()).
			product(new Object[] {
					"stream().forEach(V->out.println(V.getContig()));"
			}).
			build();
	}

	@Test(dataProvider="src2")
	public void test02(final String inputFile,final String script) 
		throws IOException
		{
		simpleTest(inputFile,script);
		}

	@DataProvider(name = "src3")
	public Object[][] createData3() {
		return new ParamCombiner().
			initList(collectAllFasta()).
			product(new Object[] {
					"stream().forEach(V->out.println(V.length()));",
					"stream().forEach(V->V.writeFasta(out));"
			}).
			build();
	}

	@Test(dataProvider="src3")
	public void test03(final String inputFile,final String script) 
		throws IOException
		{
		simpleTest(inputFile,script);
		}
	
}

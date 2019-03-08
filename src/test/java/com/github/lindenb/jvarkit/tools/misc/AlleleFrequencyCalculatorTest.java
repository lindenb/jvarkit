package com.github.lindenb.jvarkit.tools.misc;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;

public class AlleleFrequencyCalculatorTest {
	
	private final TestSupport support = new TestSupport();

	@DataProvider(name = "src1")
	public Object[][] createData1() {
		return support.toArrayArray(
				support.allVcfOrBcf().
				map(F->new Object[] {F})
				);
		}
	
	@Test(dataProvider="src1")
	public void test01(final String vcffile) throws IOException {
		try {
			final Path output = support.createTmpPath(".tsv");
			Assert.assertEquals(new AlleleFrequencyCalculator().instanceMain(new String[]{
		    		"-o",output.toString(),
		    		vcffile
		    		}),0);
			support.assertTsvTableIsConsitent(output, null);
			} finally
				{
				support.removeTmpFiles();
				}
		}
}

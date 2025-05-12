package com.github.lindenb.jvarkit.tools.misc;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;


import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.vcf.VCFUtilsTest;


public class DownSampleVcfTest {
	
	private final TestSupport support = new TestSupport();

	@DataProvider(name = "src1")
	public Object[][] createData1() {
		return support.toArrayArray(support.
				allVcfOrBcf().
				map(F->new Object[] {F})
				)
				;
		}
	
	@Test(dataProvider="src1")
    public void test01(final String vcfin) throws IOException {
		try {
	    	final Path output = support.createTmpPath(".vcf");
	    	
	        Assert.assertEquals(new DownSampleVcf().instanceMain(new String[]{
	        		"-o",output.toString(),
	        		"-n","5",
	        		vcfin
	        	}),0);
	        Assert.assertTrue(support.variantStream(output).count()<=10L);
			}
		finally
			{
			support.removeTmpFiles();
			}
    	}
}

package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.IOException;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

public class AlleleFrequencyCalculatorTest extends TestUtils{
	@Test(dataProvider="all-vcf-files")
	public void test01(final String vcffile) throws IOException {
		
		final File output = createTmpFile(".tsv");
		Assert.assertEquals(new AlleleFrequencyCalculator().instanceMain(new String[]{
	    		"-o",output.getPath(),
	    		vcffile
	    		}),0);
		super.assertTsvTableIsConsitent(output, null);
		}
}

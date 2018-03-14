package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.IOException;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;


public class DownSampleVcfTest extends TestUtils {
	@Test(dataProvider="all-vcf-files")
    public void test01(final String vcfin) throws IOException {
    	final File output = super.createTmpFile(".vcf");
    	
        Assert.assertEquals(new DownSampleVcf().instanceMain(new String[]{
        		"-o",output.getPath(),
        		"-n","5",
        		vcfin
        	}),0);
        Assert.assertTrue(variantStream(output).count()<=10L);
    	}
}

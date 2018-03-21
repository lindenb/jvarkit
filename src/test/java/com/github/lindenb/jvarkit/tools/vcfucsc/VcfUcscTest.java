package com.github.lindenb.jvarkit.tools.vcfucsc;

import java.io.File;
import java.io.IOException;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

public class VcfUcscTest extends TestUtils {
	@DataProvider(name="src01")
	public Object[][] getData() {
		return new Object[][] {
			{"vistaEnhancers","${chromStart}|${chromEnd}|${name}|${score}"},
			{"wgEncodeRegDnaseClusteredV3","${score}"},
			};
	}
	
	@Test(dataProvider="src01")
	public void test01(final String table,final String query) 
		throws IOException
		{
		final File output = super.createTmpFile(".vcf");
		
		Assert.assertEquals(0,new VcfUcsc().instanceMain(
        		newCmd().add(
        		"-o",output,
        		"--table",table,
        		"-e",query,
        		SRC_TEST_RESOURCE+"/ExAC.r1.sites.vep.vcf.gz"
        		).make()
        	));
        assertIsVcf(output);
		}

}

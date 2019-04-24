package com.github.lindenb.jvarkit.tools.vcfucsc;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;

@AlsoTest(LauncherTest.class)
public class VcfUcscTest {
	
	private final TestSupport support =new TestSupport();

	
	@DataProvider(name="src01")
	public Object[][] getData() {
		return new Object[][] {
			{"vistaEnhancers","chromStart + \"|\" + chromEnd + \"|\" + name + \"|\"+ score"},
			{"wgEncodeRegDnaseClusteredV3","score"},
			};
	}
	
	@Test(dataProvider="src01")
	public void test01(final String table,final String query) 
		throws IOException
		{
		try {
		final Path output = support.createTmpPath(".vcf");
			
			Assert.assertEquals(new VcfUcsc().instanceMain(new String[] {
	        		"-o",output.toString(),
	        		"--table",table,
	        		"-e",query,
	        		support.resource("ExAC.r1.sites.vep.vcf.gz")
					}),0);
			support.assertIsVcf(output);
			}
		finally {
			support.removeTmpFiles();
			}
		}

}

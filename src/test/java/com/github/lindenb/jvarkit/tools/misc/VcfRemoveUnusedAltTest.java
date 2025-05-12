package com.github.lindenb.jvarkit.tools.misc;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;


import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.vcf.VCFUtilsTest;


public class VcfRemoveUnusedAltTest {
	
	private final TestSupport support = new TestSupport();

	@DataProvider(name="src01")
	public Object[][] testData01() {
			return support.toArrayArray(
					support.allVcfOrBcf().
					map(S->new Object[] {S})
					);

			}
	
	@Test(dataProvider="src01")
	public void test01(final String inputFile) 
		throws IOException
			{		
			try {
				final Path output = support.createTmpPath(".vcf");
		        Assert.assertEquals(0,new VcfRemoveUnusedAlt().instanceMain(new String[]{
		        		"-o",output.toString(),
		        		"--onespan",
		        		inputFile
		        	}));
		        support.assertIsVcf(output);
				} 
			finally {
				support.removeTmpFiles();
				}
			}
	}

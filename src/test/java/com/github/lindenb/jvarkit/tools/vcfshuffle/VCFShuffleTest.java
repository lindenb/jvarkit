package com.github.lindenb.jvarkit.tools.vcfshuffle;

import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;


import com.github.lindenb.jvarkit.tools.tests.TestSupport;



public class VCFShuffleTest  {
	
	

	@DataProvider(name="src01")
	public Object[][] testData01() {
			final TestSupport support = new TestSupport();
			return support.toArrayArray(
					support.allVcfOrBcf().
					map(S->new Object[] {S})
					);
			}
	
	@Test(dataProvider="src01")
	public void test01(final String inputFile) 
		throws IOException
		{
		final TestSupport support = new TestSupport();
		try {
			final Path out = support.createTmpPath(".vcf");
			Assert.assertEquals(0,new VCFShuffle().instanceMain(new String[] {
				"-o",out.toString(),
				inputFile
				}));
			support.assertIsVcf(out);
			Assert.assertEquals(
				support.variantStream(Paths.get(inputFile)).count(),
				support.variantStream(out).count()
				);
			}
		finally
			{
			support.removeTmpFiles();
			}
		}
	}

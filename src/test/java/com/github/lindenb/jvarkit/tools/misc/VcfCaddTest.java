package com.github.lindenb.jvarkit.tools.misc;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;

public class VcfCaddTest {
	
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
			final Path out = support.createTmpPath(".vcf");
			Assert.assertEquals(0,new VcfCadd().instanceMain(new String[] {
				"-o",out.toString(),
				"-u","http://krishna.gs.washington.edu/download/CADD/v1.3/1000G_phase3.tsv.gz",
				inputFile
				}));
			support.assertIsVcf(out);
			}
		finally {
			support.removeTmpFiles();
			}
		}
}

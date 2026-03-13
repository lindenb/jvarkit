package com.github.lindenb.jvarkit.tools.vcfflatten;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;

public class VCFFlattenTest {

	

	@DataProvider(name="src01")
	public Object[][] testData01() {
			final TestSupport support = new TestSupport();
			return support.toArrayArray(
					support.allVcfOrBcf().
					filter(S->S.endsWith(".vcf.gz")).
					map(S->new Object[] {S})
					);
			}
	
	@Test(dataProvider="src01")
	public void all(final String inputFile) 
		throws IOException
		{
		final TestSupport support = new TestSupport();
		try {
			final Path out = support.createTmpPath(".vcf");
			Assert.assertEquals(0,new VCFFlatten().instanceMain(new String[] {
				"-o",out.toString(),
				inputFile
				}));
			support.assertIsVcf(out);
			}
		finally
			{
			support.removeTmpFiles();
			}
		}
	
	@Test(dataProvider="src01")
	public void sliding(final String inputFile) 
		throws IOException
		{
		final TestSupport support = new TestSupport();
		try {
			final Path out = support.createTmpPath(".vcf");
			Assert.assertEquals(0,new VCFFlatten().instanceMain(new String[] {
					"-X","+100",
					"-o",out.toString(),
				inputFile
				}));
			support.assertIsVcf(out);
			}
		finally
			{
			support.removeTmpFiles();
			}
		}
	
}

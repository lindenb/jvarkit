package com.github.lindenb.jvarkit.tools.vcf2xml;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;

public class Vcf2XmlTest {


	@DataProvider(name="src01")
	public Object[][] testData01() {
		final TestSupport support = new TestSupport();
			return support.toArrayArray(
					support.allVcfOrBcf().
					map(S->new Object[] {S})
					);
			}

	@Test(dataProvider="src01")
	public void testXml(final String inputFile) 
		throws IOException
		{
		final TestSupport support = new TestSupport();
		try {
			final Path out = support.createTmpPath(".vcf");
			final Vcf2Xml cmd =new Vcf2Xml();
			Assert.assertEquals(0,cmd.instanceMain(new String[] {
				"-o",out.toString(),
				inputFile
				}));
			support.assertIsXml(out);
			}
		finally {
			support.removeTmpFiles();
			}
		}
	@Test(dataProvider="src01")
	public void testHide(final String inputFile) 
		throws IOException
		{
		final TestSupport support = new TestSupport();
		try {
			final Path out = support.createTmpPath(".vcf");
			final Vcf2Xml cmd =new Vcf2Xml();
			Assert.assertEquals(0,cmd.instanceMain(new String[] {
				"--hide","DP,GQ,NO_CALL",
				"-o",out.toString(),
				inputFile
				}));
			support.assertIsXml(out);
			}
		finally {
			support.removeTmpFiles();
			}
		}
}

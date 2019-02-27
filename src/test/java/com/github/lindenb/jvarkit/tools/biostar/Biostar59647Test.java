package com.github.lindenb.jvarkit.tools.biostar;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;


public class Biostar59647Test {
	private final TestSupport support = new TestSupport();

	@DataProvider(name = "src1")
	public Object[][] createData1() {
		return new Object[][]{
			{support.resource("toy.bam"),support.resource("toy.fa")}
			};
	}
		
	@Test(dataProvider="src1")
	public void test1(final String bam,final String ref) throws IOException {
		try {
		final Path out = support.createTmpPath(".xml");
		Assert.assertEquals(
			new Biostar59647().instanceMain(new String[] {
			"-o",out.toString(),
			"-R",ref.toString(),
			bam}),0);
		support.assertIsXml(out);
		} 
		finally {
			support.removeTmpFiles();
			}
		}
}

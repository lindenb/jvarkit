package com.github.lindenb.jvarkit.tools.biostar;

import java.io.File;
import java.io.IOException;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;


public class Biostar59647Test extends TestUtils{
	@DataProvider(name = "src1")
	public Object[][] createData1() {
		return new Object[][]{
			{SRC_TEST_RESOURCE+"/toy.bam",SRC_TEST_RESOURCE+"/toy.fa"}
			};
	}
		
	@Test(dataProvider="src1")
	public void test1(final String bam,final String ref) throws IOException {
		final File out = createTmpFile(".xml");
		Assert.assertEquals(
			new Biostar59647().instanceMain(newCmd().
			add("-o").add(out).
			add("-R",ref,bam).
			make()
			),0);
		assertIsXml(out);
		}
}

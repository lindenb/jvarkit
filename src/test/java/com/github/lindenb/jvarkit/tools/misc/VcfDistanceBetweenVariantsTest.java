package com.github.lindenb.jvarkit.tools.misc;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;

@AlsoTest(LauncherTest.class)
public class VcfDistanceBetweenVariantsTest {
	private final TestSupport support = new TestSupport();

	@DataProvider(name = "src1")
	public Object[][] createData1() {
		return new Object[][]{
			{support.resource("toy.vcf.gz")},
			{support.resource("S1.vcf.gz")}
			};
	}
	
	@Test(dataProvider="src1")
	public void test01(final String inputFile) 
		throws IOException
		{
		try {
			final Path out = support.createTmpPath(".vcf");
			Assert.assertEquals(new VcfDistanceBetweenVariants().instanceMain(new String[] {
				"-o",out.toString(),
				inputFile
				}),0);
			support.assertIsVcf(out);
		} finally {
			support.removeTmpFiles();
			}
		}


	
}

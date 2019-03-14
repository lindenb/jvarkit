package com.github.lindenb.jvarkit.tools.biostar;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;

@AlsoTest(LauncherTest.class)
public class Biostar77288Test {
	private final TestSupport support = new TestSupport();

	
	@DataProvider(name = "src1")
	public Object[][] createData1() {
		return new Object[][]{
			{support.resource("clustalw01.msa")}
			};
	}
		
	@Test(dataProvider="src1")
	public void test01(final String url) throws IOException {
		try {
			final Path out = support.createTmpPath(".svg");
			Assert.assertEquals(
				new Biostar77288().instanceMain(new String[] {
				"-o",out.toString(), url
				}),0);
			support.assertIsXml(out);
			} 
		finally {
			support.removeTmpFiles();
			}
		}
}

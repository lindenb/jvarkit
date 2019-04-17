package com.github.lindenb.jvarkit.tools.blast;

import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;

public class BlastNToSnpTest {
	
	private final TestSupport support = new TestSupport();

	
	@DataProvider(name = "src1")
	public Object[][] createData2() {
		return new Object[][]{
			{support.resource("rotavirus_rf.blastn.01.xml")}
		};
	}

	@Test(dataProvider="src1")
	public void test01(final String blastInput) throws IOException {
		try {
			support.assertIsXml(Paths.get(blastInput));
			final Path out = support.createTmpPath(".txt");
			Assert.assertEquals(0,new BlastNToSnp().instanceMain(new String[] {
				"-o",out.toString(),
				blastInput
				}));
			}
		finally 
			{
			support.removeTmpFiles();
			}
		}


}

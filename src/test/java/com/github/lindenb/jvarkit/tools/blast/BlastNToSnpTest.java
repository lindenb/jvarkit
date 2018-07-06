package com.github.lindenb.jvarkit.tools.blast;

import java.io.File;
import java.io.IOException;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

public class BlastNToSnpTest extends TestUtils {
	@DataProvider(name = "src1")
	public Object[][] createData2() {
		return new Object[][]{
			{SRC_TEST_RESOURCE+"/rotavirus_rf.blastn.01.xml"}
		};
	}

	@Test(dataProvider="src1")
	public void test01(final String blastInput) throws IOException {
	assertIsXml(new File(blastInput));
	final File out = super.createTmpFile(".txt");
	Assert.assertEquals(0,new BlastNToSnp().instanceMain(new String[] {
		"-o",out.getPath(),
		blastInput
		}));
	
	}


}

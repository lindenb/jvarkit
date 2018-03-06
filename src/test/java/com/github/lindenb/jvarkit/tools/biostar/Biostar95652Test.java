package com.github.lindenb.jvarkit.tools.biostar;

import java.io.File;
import java.io.IOException;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

public class Biostar95652Test extends TestUtils{
	@DataProvider(name = "src1")
	public Object[][] createData1() {
		return new Object[][]{
			{"NP_077719.2  XP_513697.3 XP_001114248.1 XP_540266.3 XP_002686160.2 NP_035058.2 NP_077334.1  NP_001238962.1 NP_001108566.1"}
			};
		}
	
	@Test(dataProvider="src1")
	public void test1(final String acns) throws IOException {	
		final File out = createTmpFile(".svg");
		Assert.assertEquals(
			new Biostar95652().instanceMain(newCmd().
			add("-o").add(out).
			split(acns).
			make()
			),0);
		assertIsXml(out);
		}

}

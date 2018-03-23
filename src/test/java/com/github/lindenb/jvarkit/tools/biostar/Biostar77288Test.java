package com.github.lindenb.jvarkit.tools.biostar;

import java.io.File;
import java.io.IOException;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

public class Biostar77288Test extends TestUtils {
	@DataProvider(name = "src1")
	public Object[][] createData1() {
		return new Object[][]{
			{"http://www.tcoffee.org/Courses/Exercises/saragosa_pb_2010/practicals/practical_2/ex.1.19/file/clustalw.msa"}
			};
	}
		
	@Test(dataProvider="src1")
	public void test01(final String url) throws IOException {
		final File out = createTmpFile(".svg");
		Assert.assertEquals(
			new Biostar77288().instanceMain(newCmd().
			add("-o").add(out).add(url).make()
			),0);
		assertIsXml(out);
		}
}

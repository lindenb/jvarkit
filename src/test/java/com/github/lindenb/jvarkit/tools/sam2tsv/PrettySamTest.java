package com.github.lindenb.jvarkit.tools.sam2tsv;

import java.io.File;
import java.io.IOException;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

public class PrettySamTest extends TestUtils {
	
	@DataProvider(name = "src1")
	public Object[][] createData1() {
		return new Object[][] {
			{"./src/test/resources/toy.bam","./src/test/resources/toy.fa"}
		};
	}
	
	@Test(dataProvider="src1")
	public void test01(final String inBam,String inFasta) 
		throws IOException
		{
		final File out = super.createTmpFile(".txt");
		final PrettySam cmd =new PrettySam();
		Assert.assertEquals(cmd.instanceMain(new String[] {
			"-R",inFasta,
			"-o",out.getPath(),
			inBam
			}),0);
		assertIsNotEmpty(out);
		}
}

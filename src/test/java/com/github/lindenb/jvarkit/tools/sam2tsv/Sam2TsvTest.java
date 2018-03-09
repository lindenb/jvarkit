package com.github.lindenb.jvarkit.tools.sam2tsv;

import java.io.File;
import java.io.IOException;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

public class Sam2TsvTest extends TestUtils {
	
	@DataProvider(name = "src1")
	public Object[][] createData1() {
		return new Object[][] {
			{SRC_TEST_RESOURCE+"/toy.bam",SRC_TEST_RESOURCE+"/toy.fa"},
			{SRC_TEST_RESOURCE+"/S1.bam",SRC_TEST_RESOURCE+"/rotavirus_rf.fa"}
		};
	}
	
	
	@Test(dataProvider="src1")
	public void test01(final String inBam,String inFasta) 
		throws IOException
		{
		final File out = createTmpFile(".tsv");
		Assert.assertEquals(new Sam2Tsv().instanceMain(new String[] {
			"-R",inFasta,
			"-o",out.getPath(),
			inBam
			}),0);
		assertTsvTableIsConsitent(out,null);
		}
}

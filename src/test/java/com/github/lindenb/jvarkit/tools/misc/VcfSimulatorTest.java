package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.IOException;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

public class VcfSimulatorTest extends TestUtils{
	
	@DataProvider(name = "src1")
	public Object[][] createData1() {
		return new Object[][] {
			{SRC_TEST_RESOURCE+"/toy.fa"},
			{SRC_TEST_RESOURCE+"/rotavirus_rf.fa"}
			};
	}

@Test(dataProvider="src1")	
public void testRef(final String ref) throws IOException{
	final File out = super.createTmpFile(".vcf");
	Assert.assertEquals(0,new VcfSimulator().instanceMain(new String[] {
		"-R",ref,
		"-o",out.getPath(),
		}));
	assertIsVcf(out);
	}
}

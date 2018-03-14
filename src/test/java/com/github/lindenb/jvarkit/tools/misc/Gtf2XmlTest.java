package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.IOException;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

public class Gtf2XmlTest  extends TestUtils {
	
	@DataProvider(name="gff-data")
	public Object[][] getGffData() {
		return new Object[][] {
			{SRC_TEST_RESOURCE+"/gencode.v19.annotation.gff3","gff3"},
			{SRC_TEST_RESOURCE+"/gencode.v19.annotation.gtf","gtf"}
		};
	}
	
	@Test(dataProvider="gff-data")
	public void test01(final String inputFile,final String type) 
		throws IOException
		{
		final File out = super.createTmpFile(".xml");
		Assert.assertEquals(new Gtf2Xml().instanceMain(new String[] {
			"-o",out.getPath(),
			"--input-gtf-format",type,
			inputFile
			}),0);
		assertIsXml(out);
		}
	}

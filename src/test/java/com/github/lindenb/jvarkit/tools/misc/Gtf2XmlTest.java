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
public class Gtf2XmlTest {
	private  final TestSupport support = new TestSupport();

	
	@DataProvider(name="gff-data")
	public Object[][] getGffData() {
		return new Object[][] {
			{support.resource("gencode.v19.annotation.gff3"),"gff3"},
			{support.resource("gencode.v19.annotation.gtf"),"gtf"}
		};
	}
	
	@Test(dataProvider="gff-data")
	public void test01(final String inputFile,final String type) 
		throws IOException
		{
		try {
			final Path out = support.createTmpPath(".xml");
			Assert.assertEquals(new Gtf2Xml().instanceMain(new String[] {
				"-o",out.toString(),
				"--input-gtf-format",type,
				inputFile
				}),0);
			support.assertIsXml(out);
			}
		finally {
			support.removeTmpFiles();
			}
		}
	}

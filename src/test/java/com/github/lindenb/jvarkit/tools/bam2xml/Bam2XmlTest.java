package com.github.lindenb.jvarkit.tools.bam2xml;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;


public class Bam2XmlTest {
	
private final TestSupport support = new TestSupport();
	
	
@DataProvider(name = "src1")
public Object[][] createData1() {
	return ( Object[][])support.allSamOrBams().
			map(S->new Object[] {S}).
			toArray();
	}

@Test(dataProvider="src1")
public void test1(final String inBam) throws IOException {
	try {
		final Path out = support.createTmpPath(".xml");
		final Bam2Xml cmd =new Bam2Xml();
		Assert.assertEquals(cmd.instanceMain(new String[] {
			"-o",out.toString(),
			inBam
			}),0);
		support.assertIsXml(out);
		}
	finally
		{
		support.removeTmpFiles();
		}
	}
}

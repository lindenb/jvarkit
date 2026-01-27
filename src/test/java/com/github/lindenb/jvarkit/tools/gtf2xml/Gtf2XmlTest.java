package com.github.lindenb.jvarkit.tools.gtf2xml;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;


public class Gtf2XmlTest {
	
	@Test
	public void runSimple() throws IOException {
		Path tmp=null;
		final TestSupport support = new TestSupport();
		try {
			tmp = support.createTmpPath(".xml");
			Assert.assertEquals(
				new Gtf2Xml().instanceMain(new String[] {
					"-o",tmp.toString(),
					support.resource("Homo_sapiens.GRCh37.87.gtf.gz")
					}),0);
			support.assertIsXml(tmp);
			}
		finally {
			support.removeTmpFiles();
			}
	}
	
	
	@Test
	public void withDict() throws IOException {
		Path tmp=null;
		final TestSupport support = new TestSupport();
		try {
			tmp = support.createTmpPath(".xml");
			Assert.assertEquals(
				new Gtf2Xml().instanceMain(new String[] {
					"-R",support.resource("human_b37.dict"),
					"-o",tmp.toString(),
					support.resource("Homo_sapiens.GRCh37.87.gtf.gz")
					}),0);
			support.assertIsXml(tmp);
			}
		finally {
			support.removeTmpFiles();
			}
	}
}

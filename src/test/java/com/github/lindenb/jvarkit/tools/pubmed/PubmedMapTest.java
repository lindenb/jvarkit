package com.github.lindenb.jvarkit.tools.pubmed;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;



public class PubmedMapTest extends PubmedDumpTest {
	@Test
	public void test01() throws IOException {
		final TestSupport support = new TestSupport();
		try {
		final Path inXml = dumpAsXml(support,"Nantes lindenbaum");
		final Path out = support.createTmpPath(".xml");
		support.assertIsXml(inXml);
		Assert.assertEquals(
				new PubmedMap().instanceMain(new String[] {
				"-o",out.toString(),
				inXml.toString()
				}),0);
		support.assertIsXml(out);
		Assert.assertTrue(Files.lines(out).anyMatch(
				S->S.contains("domain=\"fr\"")));
		} finally {
			support.removeTmpFiles();
		}
	}
}

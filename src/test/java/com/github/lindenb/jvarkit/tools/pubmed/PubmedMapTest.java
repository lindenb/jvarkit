package com.github.lindenb.jvarkit.tools.pubmed;

import java.io.File;
import java.io.IOException;

import org.testng.Assert;
import org.testng.annotations.Test;

import htsjdk.samtools.util.IOUtil;

public class PubmedMapTest extends PubmedDumpTest {
	@Test
	public void test01() throws IOException {
		final File inXml = dumpAsXml("Nantes lindenbaum");
		final File out = super.createTmpFile(".xml");
		assertIsXml(inXml);
		Assert.assertEquals(
				new PubmedMap().instanceMain(newCmd().
				add("-o").add(out).
				add(inXml).
				make()
				),0);
		assertIsXml(out);
		Assert.assertTrue(IOUtil.slurpLines(out).stream().anyMatch(
				S->S.contains("domain=\"fr\"")));
		}
}

package com.github.lindenb.jvarkit.tools.pubmed;

import java.io.File;
import java.io.IOException;

import org.testng.Assert;
import org.testng.annotations.Test;

import htsjdk.samtools.util.IOUtil;


public class PubmedCodingLanguagesTest extends PubmedDumpTest {
	@Test
	public void test01() throws IOException {
		final File inXml = dumpAsXml("Python Bioinformatics 2017 Java");
		final File out = super.createTmpFile(".tsv");
		assertIsXml(inXml);
		Assert.assertEquals(
				new PubmedCodingLanguages().instanceMain(newCmd().
				add("-o").add(out).
				add(inXml).
				make()
				),0);
		Assert.assertTrue(IOUtil.slurpLines(out).stream().anyMatch(
				S->S.contains("\tpython\t")));
		}

}

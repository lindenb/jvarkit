package com.github.lindenb.jvarkit.tools.pubmed;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

import org.testng.Assert;
import org.testng.annotations.Test;

import htsjdk.samtools.util.IOUtil;

public class PubmedGenderTest extends PubmedDumpTest {
	@Test
	public void test01() throws IOException {
		final File inXml = dumpAsXml("Nantes lindenbaum");
		final File out = super.createTmpFile(".xml");
		
		final File dbFile = super.createTmpFile(".tsv");
		final PrintWriter pw = new PrintWriter(dbFile);
		pw.println("Pierre,M,1");
		pw.println("Richard,M,1");
		pw.println("Solena,F,1");
		pw.flush();
		pw.close();
		
		
		assertIsXml(inXml);
		Assert.assertEquals(
				new PubmedGender().instanceMain(newCmd().
				add("-d").add(dbFile).
				add("-o").add(out).
				add(inXml).
				make()
				),0);
		assertIsXml(out);
		Assert.assertTrue(IOUtil.slurpLines(out).stream().anyMatch(
				S->S.contains(" male=\"1")));
		}
}

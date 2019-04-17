package com.github.lindenb.jvarkit.tools.pubmed;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.AlsoTest;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;

@AlsoTest(LauncherTest.class)
public class PubmedCodingLanguagesTest extends PubmedDumpTest {
	
	@Test
	public void test01() throws IOException {
		try {
		final Path inXml = super.dumpAsXml("Python Bioinformatics 2017 Java");
		final Path out = support.createTmpPath(".tsv");
		support.assertIsXml(inXml);
		Assert.assertEquals(
				new PubmedCodingLanguages().instanceMain(new String[] {
				"-o",out.toString(),
				inXml.toString()
				}),0);
		Assert.assertTrue(Files.lines(out).anyMatch(
				S->S.contains("\tpython\t")));
		} finally {
			support.removeTmpFiles();
		}
	}

}

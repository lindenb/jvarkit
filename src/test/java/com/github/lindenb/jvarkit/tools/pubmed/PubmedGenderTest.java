package com.github.lindenb.jvarkit.tools.pubmed;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.AlsoTest;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;

@AlsoTest(LauncherTest.class)
public class PubmedGenderTest extends PubmedDumpTest {
	@Test
	public void test01() throws IOException {
		try {
			final Path inXml = super.dumpAsXml("Nantes lindenbaum");
			final Path out = support.createTmpPath(".xml");
			
			final Path dbFile = support.createTmpPath(".tsv");
			final PrintWriter pw = new PrintWriter(Files.newBufferedWriter(dbFile));
			pw.println("Pierre,M,1");
			pw.println("Richard,M,1");
			pw.println("Solena,F,1");
			pw.flush();
			pw.close();
			
			
			support.assertIsXml(inXml);
			Assert.assertEquals(
					new PubmedGender().instanceMain(new String[] {
					"-d",dbFile.toString(),
					"-o",out.toString(),
					inXml.toString()
					}),0);
			support.assertIsXml(out);
			Assert.assertTrue(Files.lines(out).anyMatch(
					S->S.contains(" male=\"1")));
		} finally 
			{
			support.removeTmpFiles();
			}
		}
}

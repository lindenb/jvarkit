package com.github.lindenb.jvarkit.tools.biostar;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;

public class Biostar81455Test {
private final TestSupport support = new TestSupport();

@Test
public void test01() throws IOException {
	try {
		final Path out = support.createTmpPath(".tsv");
		final Path in = support.createTmpPath(".txt");
		final PrintWriter pw = new PrintWriter(Files.newOutputStream(in));
		pw.println("chr22\t41258261\nchr22\t52000000\nchr22\t0");
		pw.flush();
		pw.close();
		
		Assert.assertEquals(
			new Biostar81455().instanceMain(new String[] {
			"-o",out.toString(),
			"-KG",support.resource("test_vcf01.knownGenes.txt.gz"),
			in.toString()
			}),0);
		
		support.assertTsvTableIsConsitent(out, null);
		}
	finally {
		support.removeTmpFiles();
		}
	}

}

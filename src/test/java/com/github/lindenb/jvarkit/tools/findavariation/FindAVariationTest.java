package com.github.lindenb.jvarkit.tools.findavariation;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;


public class FindAVariationTest {	

@Test
public void test01() throws IOException {
	final TestSupport support = new TestSupport();
	try {
		final Path input = support.createTmpPath(".list");
			try(PrintWriter pw=new PrintWriter(Files.newBufferedWriter(input))) {
			pw.println(support.resource("S1.vcf.gz"));
			pw.println(support.resource("S3.vcf.gz"));
			pw.println(support.resource("S4.vcf.gz"));
			pw.println(support.resource("toy.vcf.gz"));
			pw.println(support.resource("toy.bcf"));
			pw.flush();
			}
		Assert.assertTrue(support.wc(input)>0L);
		final Path output = support.createTmpPath(".tsv");
		Assert.assertEquals(new FindAVariation().instanceMain(new String[]{
	    		"-o",output.toString(),
	    		"--vcf",support.resource("S2.vcf.gz"),
	    		input.toString()
	    		}),0);
		Assert.assertTrue(support.wc(output)>1L);
		support.assertTsvTableIsConsitent(output, null);
		} finally {
			support.removeTmpFiles();
		}
	}
}

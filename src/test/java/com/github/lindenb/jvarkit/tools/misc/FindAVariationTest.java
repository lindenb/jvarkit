package com.github.lindenb.jvarkit.tools.misc;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.vcf.VCFUtilsTest;

@AlsoTest(VCFUtilsTest.class)
public class FindAVariationTest {	
	private final TestSupport support = new TestSupport();

@Test
public void test01() throws IOException {
	try {
		final Path input = support.createTmpPath(".list");
		PrintWriter pw=new PrintWriter(Files.newBufferedWriter(input));
		pw.println(support.resource("S1.vcf.gz"));
		pw.println(support.resource("S3.vcf.gz"));
		pw.println(support.resource("S4.vcf.gz"));
		pw.println(support.resource("toy.vcf.gz"));
		pw.flush();
		pw.close();
		Assert.assertTrue(support.wc(input)>0L);
		final Path output = support.createTmpPath(".tsv");
		Assert.assertEquals(new FindAVariation().instanceMain(new String[]{
	    		"-o",output.toString(),
	    		"-p","ref2:14",
	    		input.toString()
	    		}),0);
		Assert.assertTrue(support.wc(output)>1L);
		support.assertTsvTableIsConsitent(output, null);
		} finally {
			support.removeTmpFiles();
		}
	}
}

package com.github.lindenb.jvarkit.tools.minigenotyper;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;

public class MiniGenotyperTest {

	@Test
	public void testGenotype() throws IOException {
		final TestSupport support=new TestSupport();
		try {
			Path bams = support.createTmpPath(".list");
			try(PrintWriter pw  = new PrintWriter(Files.newBufferedWriter(bams))) {
				pw.println(support.resource("S1.bam"));
				pw.println(support.resource("S2.bam"));
				pw.println(support.resource("S3.bam"));
				pw.println(support.resource("S4.bam"));
				pw.flush();
				}
			Path out = support.createTmpPath(".vcf");

			Assert.assertEquals(new MiniGenotyper().instanceMain(new String[] {
					"-o",out.toString(),
					"-V",support.resource("rotavirus_rf.vcf.gz"),
					"-R",support.resource("rotavirus_rf.fa"),
					bams.toString()
					}),0);
			support.assertIsVcf(out);
			}
		finally {
			support.removeTmpFiles();
		}
	}

}

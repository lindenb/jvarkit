package com.github.lindenb.jvarkit.tools.structvar;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.io.IOUtils;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;


public class ValidateCnvTest {
	final TestSupport support = new TestSupport();
	@Test
	public void test01() throws IOException
		{
		try {
		final Path bamlist= support.createTmpPath(".list");
		final PrintWriter pw = IOUtils.openPathForPrintWriter(bamlist);
		pw.println(support.resource("S1.bam"));
		pw.println(support.resource("S2.bam"));
		pw.println(support.resource("S3.bam"));
		pw.println(support.resource("S4.bam"));
		pw.println(support.resource("S5.bam"));
		pw.flush();
		pw.close();
		
		
		final Path out = support.createTmpPath(".vcf");
		Assert.assertEquals(new ValidateCnv().instanceMain(new String[] {
				"-o",out.toString(),
				"-B",bamlist.toString(),
				"-R",support.resource("rotavirus_rf.fa"),
			support.resource("rotavirus_rf.vcf.gz")
			}),0
			);
		support.assertIsVcf(out);
		} finally {
			support.removeTmpFiles();
			}
		}
	}

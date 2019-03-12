package com.github.lindenb.jvarkit.tools.structvar;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.tools.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;

@AlsoTest(LauncherTest.class)
public class ValidateCnvTest {
	final TestSupport support = new TestSupport();
	@Test
	public void test01() throws IOException
		{
		try {
		final Path bed= support.createTmpPath(".bed");
		final PrintWriter pw = IOUtils.openPathForPrintWriter(bed);
		for(int i=0;i< 10;++i) {
			int p = support.random.nextInt(100);
			int L = support.random.nextInt(100);
			pw.println("RF01"+i+"\t"+p+"\t"+(p+L)+"\tid"+i);
			}
		pw.flush();
		pw.close();
		
		
		final Path out = support.createTmpPath(".vcf");
		Assert.assertEquals(new ValidateCnv().instanceMain(new String[] {
				"-o",out.toString(),
				"-B",bed.toString(),
			support.resource("S1.bam"),
			support.resource("S2.bam"),
			support.resource("S3.bam"),
			support.resource("S4.bam")
			}),0
			);
		support.assertIsVcf(out);
		} finally
		{
			support.removeTmpFiles();
		}
		}
	}

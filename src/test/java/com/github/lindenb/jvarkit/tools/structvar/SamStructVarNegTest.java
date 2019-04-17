package com.github.lindenb.jvarkit.tools.structvar;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;

@AlsoTest(LauncherTest.class)
public class SamStructVarNegTest {
	final TestSupport support = new TestSupport();

	@Test
	public void test01() 
		throws IOException
		{
		try {
			final Path list = support.createTmpPath(".list");
			PrintWriter pw = new PrintWriter(Files.newBufferedWriter(list));
			pw.println(support.resource("S2.bam"));
			pw.println(support.resource("S3.bam"));
			pw.println(support.resource("S4.bam"));
			pw.println(support.resource("S5.bam"));
			pw.flush();
			pw.close();
			
			final Path out = support.createTmpPath(".vcf");
			Assert.assertEquals(new SamStructVarNeg().instanceMain(new String[] {
					"-o",out.toString(),
					"--bams",list.toString(),
				support.resource("S1.bam")
				}),0);
			support.assertIsVcf(out);
		} finally {
			support.removeTmpFiles();
			}
		}
}

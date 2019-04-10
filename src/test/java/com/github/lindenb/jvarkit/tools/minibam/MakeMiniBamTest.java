package com.github.lindenb.jvarkit.tools.minibam;

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
public class MakeMiniBamTest {
	final TestSupport support= new TestSupport();

	@Test
	public void test01() throws IOException {
		try {
			Path bams = support.createTmpPath(".list");
			PrintWriter pw  = new PrintWriter(Files.newBufferedWriter(bams));
			pw.println(support.resource("S1.bam"));
			pw.println(support.resource("S2.bam"));
			pw.println(support.resource("S3.bam"));
			pw.println(support.resource("S4.bam"));
			pw.flush();
			pw.close();
			Path out = support.createTmpPath(".zip");

			Assert.assertEquals(new MakeMiniBam().instanceMain(new String[] {
					"-o",out.toString(),
					"-p","RF01:100",
					"-p","RF02:100",
					bams.toString()
					}),0);
			support.assertZip(out);
			}
		finally {
			support.removeTmpFiles();
		}
	}
}

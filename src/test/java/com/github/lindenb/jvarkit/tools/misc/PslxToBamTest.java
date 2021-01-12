package com.github.lindenb.jvarkit.tools.misc;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;


@AlsoTest({LauncherTest.class})
public class PslxToBamTest {
	private final TestSupport support = new TestSupport();
	@Test
	public void test01() throws IOException {
		try {
			final Path in= support.createTmpPath(".psl");
			final String lines="101\t0\t0\t0\t1\t29\t0\t0\t-\tRF01:100-200/rc+N+INS\t140\t5\t135\tRF01\t3302\t99\t200\t2\t67,34,\t5,101,\t99,166,\ttattcttccaatagtgaattagagaatagatgtattgaatttcattctaaatgcttagaaaactcaa,agaatggactatcattgaaaaagctctttgttga,\ttattcttccaatagtgaattagagaatagatgtattgaatttcattctaaatgcttagaaaactcaa,agaatggactatcattgaaaaagctctttgttga,\n";
			Files.write(in,lines.getBytes());
			final Path out = support.createTmpPath(".bam");
			Assert.assertEquals(
				new PslxToBam().instanceMain(new String[] {
					"-o",out.toString(),
					"-R",support.resource("rotavirus_rf.fa"),
					in.toString()
				}),0);
			support.assertIsValidBam(out);
		} finally 
		{
			support.removeTmpFiles();
		}
	}
	
	
}
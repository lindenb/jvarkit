package com.github.lindenb.jvarkit.tools.fastq;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;

@AlsoTest(LauncherTest.class)
public class FastqShuffleTest {
	private final TestSupport support = new TestSupport();
	
	@Test
	public void test01() throws IOException {
		try {
			final Path out = support.createTmpPath(".fq");
			Assert.assertEquals(new FastqShuffle().instanceMain(new String[] {
				"-o",out.toString(),
				support.resource("S1.R1.fq.gz"),
				support.resource("S1.R2.fq.gz"),
				}),0);
			support.assertIsFastq(out);
		} finally {
			support.removeTmpFiles();
		}
	}

}

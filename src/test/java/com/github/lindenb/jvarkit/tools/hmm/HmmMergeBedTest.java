package com.github.lindenb.jvarkit.tools.hmm;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;

@AlsoTest(Launcher.class)
public class HmmMergeBedTest {
	private final TestSupport support = new TestSupport();

	@Test
	public void test01() throws IOException {
		Path tmp1=null;
		Path tmp2=null;
		Path tmp3=null;
		Path out=null;
		try {
			tmp1 = support.createTmpPath(".bed");
			Files.write(tmp1, Arrays.asList("chr1\t10\t100\tx1"));
			tmp2 = support.createTmpPath(".bed");
			Files.write(tmp2, Arrays.asList("chr1\t20\t110\tx1"));
			tmp3 = support.createTmpPath(".bed");
			Files.write(tmp2, Arrays.asList("chr1\t18\t110\tx3"));
			
			out = support.createTmpPath(".bed");

			
			
			Assert.assertEquals(new HmmMergeBed().instanceMain(new String[] {
					"-o",out.toString(),
					"-n","2",
					tmp1.toString(),
					tmp2.toString(),
					tmp3.toString()
					}),0);
			support.assertIsBed(out);
			}
		finally {
			support.removeTmpFiles();
			}
	}
}

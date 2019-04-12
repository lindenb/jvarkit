package com.github.lindenb.jvarkit.tools.liftover;

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
public class ConvertLiftOverChainTest {
	final TestSupport support= new TestSupport();

	@Test
	public void test01() throws IOException {
		try {
			Path chain = support.createTmpPath(".chain");
			PrintWriter pw  = new PrintWriter(Files.newBufferedWriter(chain));
			pw.println("chain 20851231461 chr1 249250621 + 10000 249240621 chr1 248956422 + 10000 248946422 2");
			pw.println("chain 2920 chrUn_gl000228 129120 + 104295 104331 chr10 133797422 + 133748242 133748278 1495");
			pw.flush();
			pw.close();
			Path out = support.createTmpPath(".chain");

			Assert.assertEquals(new ConvertLiftOverChain().instanceMain(new String[] {
					"-o",out.toString(),
					"-R1",support.resource("human_b37.dict"),
					"-R2",support.resource("human_b37.dict"),
					chain.toString()
					}),0);
			Assert.assertTrue(Files.lines(out).noneMatch(L->L.contains("chr")));
			}
		finally {
			support.removeTmpFiles();
		}
	}
}

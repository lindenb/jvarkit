package com.github.lindenb.jvarkit.tools.sam2tsv;

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
public class WesCnvTViewTest {
	private final TestSupport support = new TestSupport();

	
	@Test
	public void test01() throws IOException {
		try {
		final  Path bamList = support.createTmpPath(".list");
		final PrintWriter pw = new PrintWriter(Files.newBufferedWriter(bamList));
		for(int i=1;i<=5;i++)
			{
			pw.println(support.resource("S"+i+".bam"));
			}
		pw.flush();
		pw.close();
		final Path out = support.createTmpPath(".txt");
		Assert.assertEquals(new WesCnvTView().instanceMain(new String[] {
				"-l",bamList.toString(),
				"-o",out.toString(),
				"RF01:200-1000",
				"RF02:300-500"
				}),0);
		support.assertIsNotEmpty(out);
		}
	finally {
		support.removeTmpFiles();
		}
	}
}

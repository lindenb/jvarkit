package com.github.lindenb.jvarkit.tools.sam2tsv;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;

@AlsoTest(LauncherTest.class)
public class CnvTViewTest {
	private final TestSupport support = new TestSupport();

	
	@Test
	public void test01() throws IOException {
		try {
		final  Path bamList = support.createTmpPath(".list");
		try(PrintWriter pw = new PrintWriter(Files.newBufferedWriter(bamList))) {
			for(int i=1;i<=5;i++)
				{
				pw.println(support.resource("S"+i+".bam"));
				}
			pw.flush();
			}
		final Path out = support.createTmpPath(".txt");
		Assert.assertEquals(new CnvTView().instanceMain(new String[] {
				"-o",out.toString(),
				"-R " ,support.resource("rotavirus_rf.fa"),
				"-r " ,"RF01:200-1000",
				bamList.toString()
				}),0);
		support.assertIsNotEmpty(out);
		}
	finally {
		support.removeTmpFiles();
		}
	}
}

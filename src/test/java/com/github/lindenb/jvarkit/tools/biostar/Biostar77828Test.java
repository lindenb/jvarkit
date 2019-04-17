package com.github.lindenb.jvarkit.tools.biostar;

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
public class Biostar77828Test {
@Test
public void test01() throws IOException {
	final TestSupport support = new TestSupport();
	try {
		Path bed = support.createTmpPath(".bed");
		PrintWriter pw=new PrintWriter(Files.newBufferedWriter(bed));
		for(int i=0;i< 1000;i++)
			{
			int L= 1+support.random.nextInt(1000);
			int p= support.random.nextInt(100000);
			pw.println("chr"+(i%10)+"\t"+p+"\t"+(p+L));
			}
		pw.flush();
		pw.close();
		support.assertIsBed(bed);
		
		Path out = support.createTmpPath(".txt");

		Assert.assertEquals(
				new Biostar77828().instanceMain(new String[] {
				"-o",out.toString(),
				"--iter","1000",
				bed.toString()
				}),0);
		support.assertIsNotEmpty(out);
	} finally {
		support.removeTmpFiles();
	}
	
	}
}

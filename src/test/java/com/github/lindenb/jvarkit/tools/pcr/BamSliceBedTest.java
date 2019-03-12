package com.github.lindenb.jvarkit.tools.pcr;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;


public class BamSliceBedTest  {
	private final TestSupport support = new TestSupport();

	
	@Test
	public void test01() throws IOException {
		try {
			final Path out = support.createTmpPath(".bam");
			final Path bed = support.createTmpPath(".bed");
			final PrintWriter pw = new PrintWriter(Files.newBufferedWriter(bed));
			pw.println("RF01\t1\t10");
			pw.println("RF01\t5\t15");
			pw.println("RF02\t15\t20");
			pw.flush();
			pw.close();
			support.assertIsBed(bed);
			
			Assert.assertEquals(new BamSliceBed().instanceMain(new String[]{
		    		"-o",out.toString(),
		    		"-B",bed.toString(),
		    		support.resource("S1.bam")
		    		}),0);
			support.assertIsValidBam(out);
			}
		finally {
			support.removeTmpFiles();
			}
		}

}

package com.github.lindenb.jvarkit.tools.pcr;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;


public class BamSliceBedTest extends TestUtils {
	
	
	@Test
	public void test01() throws IOException {
		final File out = super.createTmpFile(".bam");
		final File bed = super.createTmpFile(".bed");
		final PrintWriter pw = new PrintWriter(bed);
		pw.println("RF01\t1\t10");
		pw.println("RF01\t5\t15");
		pw.println("RF02\t15\t20");
		pw.flush();
		pw.close();
		super.assertIsBed(bed);
		
		Assert.assertEquals(new BamSliceBed().instanceMain(new String[]{
	    		"-o",out.getPath(),
	    		"-V",out.getPath(),
	    		SRC_TEST_RESOURCE+"/S1.bam"
	    		}),0);
		assertIsValidBam(out);
		}

}

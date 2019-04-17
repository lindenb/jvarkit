package com.github.lindenb.jvarkit.tools.xcontamination;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;

@AlsoTest(LauncherTest.class)
public class XContaminationsTest {
	
	private final TestSupport support =new TestSupport();

@Test	
public void test01() throws IOException {
	try {
		Path output = support.createTmpPath(".vcf");
		Assert.assertEquals(new XContaminations().instanceMain(new String[] {
	    		"-ov","-sample",
	    		"-o",output.toString(),
	    		support.resource("rotavirus_rf.vcf.gz"),
	    		support.resource("S1.bam"),
	    		support.resource("S2.bam"),
	    		support.resource("S3.bam"),
	    		support.resource("S4.bam"),
	    		support.resource("S5.bam")
				}),0);
		support.assertIsVcf(output);
		} 
	finally 
		{
		support.removeTmpFiles();
		}
	}
}

package com.github.lindenb.jvarkit.tools.vcfanonymize;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;

public class VcfAnonymizeTest {
	
	@Test
	public void testAnonymize() throws IOException
		{
		final TestSupport support=new TestSupport();
		try {
			final String inputFile = support.resource("rotavirus_rf.vcf.gz");
			
		
			final Path output = support.createTmpPath(".vcf");
			
			Assert.assertEquals(new VcfAnonymize().instanceMain(new String[] {
	        		"-o",output.toString(),
	        		inputFile.toString()}),0);
			support.assertIsVcf(output);
			} 
		finally {
			support.removeTmpFiles();
			}
		}

}

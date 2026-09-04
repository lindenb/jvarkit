package com.github.lindenb.jvarkit.tools.vcfcoocurrence;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;

public class VcfCoOccurrenceTest {

	
	
	@Test
	public void testFile(final String inputFile,String expr) 
		throws IOException
		{
		final TestSupport support =new TestSupport();
		try {
			final Path output = support.createTmpPath(".tsv");
	        Assert.assertEquals(new VcfCoOccurrence().instanceMain(new String[] {
	        		"-o",output.toString(),
	        		"--cases","S1,S2",
	        		"--controls","S3,S4",
	        		support.resource("rotavirus_rf.vcf.gz")}),
	        		0);
	        support.assertTsvTableIsConsitent(output, null);
			} 
		finally
			{	
			support.removeTmpFiles();
			}
		}

}
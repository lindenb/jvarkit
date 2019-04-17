package com.github.lindenb.jvarkit.tools.vcfbigwig;

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
public class VCFBigWigTest  {
	
	private final TestSupport support =new TestSupport();

	
@Test
public void testLocal() throws IOException {
	try {
	Path out = support.createTmpPath(".vcf");
	Assert.assertEquals(0,new VCFBigWig().instanceMain(new String[] {
			"-B",support.resource("Uniqueness35bp.bigWig"),
			"-T","XXX",
			"-o",out.toString(),
			support.resource("test_vcf01.vcf")
			}));
	support.assertIsVcf(out);
	Assert.assertTrue(support.variantStream(out).anyMatch(P->P.hasAttribute("XXX")));
	} finally
	{
	support.removeTmpFiles();	
	}
}
@Test
public void testXml() throws IOException {
	try {
		Path xml = support.createTmpPath(".xml");
		PrintWriter pw = new PrintWriter(Files.newBufferedWriter(xml));
		pw.println("<registry><bigwig>");
		pw.println("<uri>"+ support.resource("Uniqueness35bp.bigWig")+"</uri>");
		pw.println("<tag>XXX</tag>");
		pw.println("<description>XXX</description>");
		pw.println("</bigwig></registry>");
		pw.flush();
		pw.close();
		support.assertIsXml(xml);
		
		Path out = support.createTmpPath(".vcf");
		Assert.assertEquals(0,new VCFBigWig().instanceMain(new String[] {
				"-B",xml.toString(),
				"-o",out.toString(),
				support.resource("test_vcf01.vcf")
				}));
		support.assertIsVcf(out);
		Assert.assertTrue(support.variantStream(out).anyMatch(P->P.hasAttribute("XXX")));
		} 
	finally
		{
		support.removeTmpFiles();
		}
	}
}

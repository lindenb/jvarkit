package com.github.lindenb.jvarkit.tools.vcfbigwig;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

public class VCFBigWigTest extends TestUtils {
@Test
public void testLocal() throws IOException {
	File out = super.createTmpFile(".vcf");
	Assert.assertEquals(0,new VCFBigWig().instanceMain(new String[] {
			"-B",SRC_TEST_RESOURCE+"/Uniqueness35bp.bigWig",
			"-T","XXX",
			"-o",out.getPath(),
			SRC_TEST_RESOURCE+"/test_vcf01.vcf"
			}));
	assertIsVcf(out);
	Assert.assertTrue(variantStream(out).anyMatch(P->P.hasAttribute("XXX")));
	}
@Test
public void testXml() throws IOException {
	File xml = super.createTmpFile(".xml");
	PrintWriter pw = new PrintWriter(xml);
	pw.println("<registry><bigwig>");
	pw.println("<uri>"+SRC_TEST_RESOURCE+"/Uniqueness35bp.bigWig</uri>");
	pw.println("<tag>XXX</tag>");
	pw.println("<description>XXX</description>");
	pw.println("</bigwig></registry>");
	pw.flush();
	pw.close();
	assertIsXml(xml);
	
	File out = super.createTmpFile(".vcf");
	Assert.assertEquals(0,new VCFBigWig().instanceMain(new String[] {
			"-B",xml.getPath(),
			"-o",out.getPath(),
			SRC_TEST_RESOURCE+"/test_vcf01.vcf"
			}));
	assertIsVcf(out);
	Assert.assertTrue(variantStream(out).anyMatch(P->P.hasAttribute("XXX")));
	}

}

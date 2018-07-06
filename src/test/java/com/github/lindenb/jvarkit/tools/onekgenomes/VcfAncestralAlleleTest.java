package com.github.lindenb.jvarkit.tools.onekgenomes;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

import org.testng.Assert;
import org.testng.annotations.Test;
import org.testng.annotations.DataProvider;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

import htsjdk.samtools.util.IOUtil;

public class VcfAncestralAlleleTest extends TestUtils {
	@DataProvider(name = "src1")
	public Object[][] createData1() {
		return new Object[][] {
			{SRC_TEST_RESOURCE+"/rotavirus_rf.fa",SRC_TEST_RESOURCE+"/rotavirus_rf.vcf.gz"},
			{SRC_TEST_RESOURCE+"/toy.fa",SRC_TEST_RESOURCE+"/toy.vcf.gz"}
		};
	}
	
	@Test(dataProvider="src1")
	public void test01(
			final String fasta,
			final String inputFile
			) throws IOException
	{
	final File output = super.createTmpFile(".vcf");
	final File manifest = super.createTmpFile(".MF");
	final File ref=new File(fasta);
	PrintWriter pw = new PrintWriter(manifest);
	for(final String l: IOUtil.slurpLines( new File(ref.getParentFile(),ref.getName()+".fai"))) {
		final String tokens[] = l.split("[\t]");
		pw.println(tokens[0]+"|"+tokens[0]+"xxx\t"+tokens[0]+"\t"+ref);
	}
	pw.flush();
	pw.close();
	 Assert.assertEquals(new VcfAncestralAllele().instanceMain(
     		newCmd().add(
     		"-o",output.getPath(),
     		"-m",manifest,
     		inputFile
     		).make()),0);
	assertIsVcf(output);
	
	Assert.assertTrue(variantStream(output).
			allMatch(V->V.hasAttribute("AA") && V.getReference().getDisplayString().equalsIgnoreCase(V.getAttributeAsString("AA", "."))));
	}
}

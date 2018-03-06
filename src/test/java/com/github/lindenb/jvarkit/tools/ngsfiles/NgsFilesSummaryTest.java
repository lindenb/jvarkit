package com.github.lindenb.jvarkit.tools.ngsfiles;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

import htsjdk.samtools.util.IOUtil;

public class NgsFilesSummaryTest extends TestUtils{
@Test
public void test01()throws IOException {
	final File pathFile = super.createTmpFile(".txt");
	final PrintWriter pw = new PrintWriter(pathFile);
	_collectFiles(new File("./src/test/resources/"),
			(D,N)->N.endsWith(".vcf") || 
			       N.endsWith(".vcf.gz") || 
			       N.endsWith(".bam")  || 
			       N.endsWith(".fq.gz")
			       ).forEach(L->pw.println(L));
	pw.flush();
	pw.close();
	final File out =  super.createTmpFile(".txt");
	Assert.assertEquals(new NgsFilesSummary().instanceMain(newCmd().
		add("-o",out.getPath()).
		add(pathFile).make()
		),0);
	Assert.assertTrue(IOUtil.slurpLines(out).stream().
			anyMatch(L->L.contains("VCF") && 
					L.contains("S1") && 
					L.contains("vcf.gz"))
			);
	}
}

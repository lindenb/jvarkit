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
	assertTsvTableIsConsitent(out, null);
	Assert.assertTrue(IOUtil.slurpLines(out).stream().
			anyMatch(L->L.contains("VCF") && 
					L.contains("S1") && 
					L.contains("vcf.gz"))
			);
	}

public void test02() 
		throws IOException
		{
		final File tmp1 = super.createTmpFile(".txt");
		final PrintWriter pw = new PrintWriter(tmp1);
		for(int i=1;i<=5;++i)
			{
			pw.println(SRC_TEST_RESOURCE+"/S"+i+".bam");
			pw.println(SRC_TEST_RESOURCE+"/S"+i+".vcf.gz");
			}
		pw.println(SRC_TEST_RESOURCE+"/rotavirus_rf.vcf.gz");
		pw.flush();
		pw.close();
		
		final File out = super.createTmpFile(".txt");
		Assert.assertEquals(0,new NgsFilesSummary().instanceMain(new String[] {
			"-o",out.getPath(),
			tmp1.getPath()
			}));
		assertTsvTableIsConsitent(out, null);
		}
}

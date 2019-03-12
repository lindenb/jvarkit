package com.github.lindenb.jvarkit.tools.ngsfiles;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;
public class NgsFilesSummaryTest {
	
private	final TestSupport support = new TestSupport();	
	
@Test
public void test01()throws IOException {
	try {
		final Path pathFile = support.createTmpPath(".txt");
		final PrintWriter pw = new PrintWriter(Files.newBufferedWriter(pathFile));
		support.allVcfOrBcf().forEach(L->pw.println(L));
		support.allSamOrBams().forEach(L->pw.println(L));
		pw.flush();
		pw.close();
		final Path out =  support.createTmpPath(".txt");
		Assert.assertEquals(new NgsFilesSummary().instanceMain(new String[] {
			"-o",out.toString(),
			pathFile.toString()
			}),0);
		support.assertTsvTableIsConsitent(out, null);
		Assert.assertTrue(Files.lines(out).
				anyMatch(L->L.contains("VCF") && 
						L.contains("S1") && 
						L.contains("vcf.gz"))
				);
		} 
	finally
		{
		support.removeTmpFiles();
		}
	}

@Test
public void test02() 
		throws IOException
		{
		try {
			final Path tmp1 = support.createTmpPath(".txt");
			final PrintWriter pw = new PrintWriter(Files.newBufferedWriter(tmp1));
			for(int i=1;i<=5;++i)
				{
				pw.println(support.resource("S"+i+".bam"));
				pw.println(support.resource("S"+i+".vcf.gz"));
				}
			pw.println(support.resource("rotavirus_rf.vcf.gz"));
			pw.flush();
			pw.close();
			
			final Path out = support.createTmpPath(".txt");
			Assert.assertEquals(0,new NgsFilesSummary().instanceMain(new String[] {
				"-o",out.toString(),
				tmp1.toString()
				}));
			support.assertTsvTableIsConsitent(out, null);
			}
		finally {
			support.removeTmpFiles();
			}
		}
}

package com.github.lindenb.jvarkit.tools.bamstats04;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Random;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;

public class BamStats04Test  {
	
private final TestSupport support = new TestSupport();
private final Random random=new Random(System.currentTimeMillis());
	
@DataProvider(name = "src1")
public Object[][] createData1() {
		return (Object[][])support.allIndexedBams().
				map(S->SAMSequenceDictionaryExtractor.extractDictionary(Paths.get(S))!=null).
				map(S->new Object[] {S}).toArray();
		}
	
@Test(dataProvider="src1")
public void test1(final String inBam) throws IOException {
	try {
		final Path out = support.createTmpPath(".txt");
		final Path bedout = support.createTmpPath(".bed");
		final SAMSequenceDictionary dict= SequenceDictionaryUtils.extractRequired(Paths.get(inBam));
		final PrintWriter pw = new PrintWriter(Files.newOutputStream(bedout));
		for(int i=0;i<100;i++)
			{	
			final SAMSequenceRecord ssr=dict.getSequence(this.random.nextInt(dict.size()));
			int start = this.random.nextInt(ssr.getSequenceLength()-1);
			int end = Math.min(ssr.getSequenceLength(),start+1+this.random.nextInt(100));
			pw.println(ssr.getSequenceName()+"\t"+start+"\t"+end);
			}
		pw.flush();
		pw.close();
		
	
		final BamStats04 cmd =new BamStats04();
		Assert.assertEquals(cmd.instanceMain(new String[] {
			"-o",out.toString(),
			"--bed",bedout.toString(),
			inBam
			}),0);
		support.assertIsNotEmpty(out);
		} 
	finally {
		support.removeTmpFiles();
		}
	}
}

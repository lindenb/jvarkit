package com.github.lindenb.jvarkit.tools.bamstats04;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashSet;
import java.util.Set;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;

public class BamStats05Test extends TestUtils {
@DataProvider(name = "src1")
public Object[][] createData1() {
	return new ParamCombiner().
		initList(collectIndexedBams()).
		build();
		}
	
@Test(dataProvider="src1")
public void test1(final String inBam) throws IOException {
	final File out = createTmpFile(".txt");
	final File bedout = createTmpFile(".bed");
	final SAMSequenceDictionary dict= SAMSequenceDictionaryExtractor.extractDictionary(new File(inBam));
	if(dict==null) return;
	final PrintWriter pw = new PrintWriter(bedout);
	final Set<String> seen_contig = new HashSet<>();
	for(int i=0;i<100;i++)
		{
		final SAMSequenceRecord ssr=dict.getSequence(this.random.nextInt(dict.size()));
		if(ssr.getSequenceLength()> 10_000_000) continue;
		if(seen_contig.contains(ssr.getSequenceName())) continue;
		seen_contig.add(ssr.getSequenceName());
		int start = this.random.nextInt(ssr.getSequenceLength()-1);
		int end = Math.min(ssr.getSequenceLength(),start+1+this.random.nextInt(100));
		pw.println(ssr.getSequenceName()+"\t"+start+"\t"+end+"\tGene"+ssr.getSequenceName());
		}
	pw.flush();
	pw.close();
	if(seen_contig.isEmpty()) return;//contigs too large for test
	
	Assert.assertEquals(new BamStats05().instanceMain(new String[] {
		"-o",out.getPath(),
		"--bed",bedout.getPath(),
		inBam
		}),0);
	assertTsvTableIsConsitent(out, null);
	}
}

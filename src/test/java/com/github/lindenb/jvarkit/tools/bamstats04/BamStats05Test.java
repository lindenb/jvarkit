package com.github.lindenb.jvarkit.tools.bamstats04;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;

@AlsoTest(LauncherTest.class)

public class BamStats05Test {
private final TestSupport support = new TestSupport();
private final Random random=new Random(System.currentTimeMillis());

@DataProvider(name = "src1")
public Object[][] createData1() {
	return support.toArrayArray(support.allIndexedBams().
			map(S->support.getReferenceRegistry().getReferenceByPath(Paths.get(S)).isPresent()).
			map(S->new Object[] {S})
			);
	}
	
@Test(dataProvider="src1")
public void test1(final String inBam) throws IOException {
	try {
	final Path out = support.createTmpPath(".txt");
	final Path bedout = support.createTmpPath(".bed");
	final SAMSequenceDictionary dict= SAMSequenceDictionaryExtractor.extractDictionary(Paths.get(inBam));
	if(dict==null) return;
	
	final PrintWriter pw = new PrintWriter(Files.newOutputStream(bedout));
	final Set<String> seen_contig = new HashSet<>();
	for(int i=0;i<100;i++)
		{
		final SAMSequenceRecord ssr=dict.getSequence(this.random.nextInt(dict.size()));
		if(ssr.getSequenceLength()< 100) continue;
		if(seen_contig.contains(ssr.getSequenceName())) continue;
		seen_contig.add(ssr.getSequenceName());
		int start = this.random.nextInt(ssr.getSequenceLength()-1);
		int end = Math.min(ssr.getSequenceLength(),start+1+this.random.nextInt(100));
		pw.println(ssr.getSequenceName()+"\t"+start+"\t"+end+"\tGene"+ssr.getSequenceName());
		}
	pw.flush();
	pw.close();
	if(seen_contig.isEmpty()) return;//contigs too large for test
	
	final List<String> args = new ArrayList<>();
	args.add("-o");
	args.add(out.toString());
	args.add("--bed");
	args.add(bedout.toString());
	
	
	Assert.assertEquals(new BamStats05().instanceMain(args),0);
	support.assertTsvTableIsConsitent(out, null);
	} finally {
		support.removeTmpFiles();
	}
	}
}

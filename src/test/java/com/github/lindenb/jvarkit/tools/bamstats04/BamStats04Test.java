package com.github.lindenb.jvarkit.tools.bamstats04;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;

public class BamStats04Test extends TestUtils {
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
		"-o",out.getPath(),
		"--bed",bedout.getPath(),
		inBam
		}),0);
	assertIsNotEmpty(out);
	}
}

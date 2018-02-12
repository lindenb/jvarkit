package com.github.lindenb.jvarkit.tools.biostar;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashSet;
import java.util.Random;
import java.util.Set;
import java.util.stream.Collectors;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.StringUtil;

public class Biostar214299Test extends TestUtils {
	@DataProvider(name = "src1")
	public Object[][] createData2() {
		return new ParamCombiner().
			initList(collectAllSamOrBam()).
			build();
	}

	@Test(dataProvider="src1")
	public void test01(final String inputBam) throws IOException {
	SamReader sr  = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.LENIENT).open(SamInputResource.of(inputBam));
	SAMFileHeader header = sr.getFileHeader();
	sr.close();
	Set<String> samples = header.getReadGroups().stream().
			map(S->S.getSample()).
			filter(S->!StringUtil.isBlank(S)).collect(Collectors.toSet());
	if(samples.isEmpty()) return;
	
	final SAMSequenceDictionary dict = header.getSequenceDictionary();
	final SAMSequenceRecord ssr = dict.getSequence(0);
	if(dict==null || dict.isEmpty()) return;
	final File table = File.createTempFile(".tmp.", ".txt");
	final Random random = new Random();
	final PrintWriter pw = new PrintWriter(table);
	final Set<Integer> positions = new HashSet<>();
	for(int i=0;i< 30;i++) {
		positions.add(1+random.nextInt(ssr.getSequenceLength()));
		}
	for(final Integer p:positions) {
		for(final String sample:samples) {
			pw.print(ssr.getSequenceName());
			pw.print('\t');
			pw.print(p);
			pw.print('\t');
	 		switch(random.nextInt(4))
	 			{
	 			case 0:	pw.print('A'); break;
	 			case 1:	pw.print('C'); break;
	 			case 2:	pw.print('G'); break;
	 			default:	pw.print('T'); break;
	 			}
	 		pw.print('\t');
			pw.print(sample);
			pw.println();
			}
	}
	pw.flush();
	pw.close();
	
	final File out = File.createTempFile(".tmp.", ".bam");
	final Biostar214299  cmd =new Biostar214299();
	Assert.assertEquals(0,cmd.instanceMain(new String[] {
		"-p",table.getPath(),
		"-o",out.getPath(),
		inputBam
		}));
	Assert.assertTrue(out.delete());
	Assert.assertTrue(table.delete());

	}
}

package com.github.lindenb.jvarkit.tools.biostar;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.HashSet;
import java.util.Random;
import java.util.Set;
import java.util.stream.Collectors;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.StringUtil;

public class Biostar214299Test {
	
	private final TestSupport support = new TestSupport();

	
	@DataProvider(name = "src1")
	public Object[][] createData1() {
		return (Object[][])support.
				allSamOrBams().
				map(F->new Object[] {F}).
				toArray()
				;
		}


	@Test(dataProvider="src1")
	public void test01(final String inputBam) throws IOException {
		try {
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
		final Path table = support.createTmpPath(".txt");
		final Random random = new Random();
		final PrintWriter pw = new PrintWriter(Files.newOutputStream(table));
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
		
		final Path out = support.createTmpPath(".bam");
		Assert.assertEquals(new Biostar214299().instanceMain(new String[] {
			"-p",table.toString(),
			"-o",out.toString(),
			inputBam
			}),0);
		support.assertIsValidBam(out);
		} 
		finally {
		support.removeTmpFiles();
		}
	}
}

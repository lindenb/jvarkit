package com.github.lindenb.jvarkit.tools.misc;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.stream.Collectors;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceDictionaryCodec;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;

public class VcfSetSequenceDictionaryTest {
	
	private final TestSupport support = new TestSupport();

	@DataProvider(name="src01")
	public Object[][] testData01() {
			return support.toArrayArray(
					support.allVcfOrBcf().
					map(S->new Object[] {S})
					);
			}
	
	private Path prepare(final String inputFile)  throws IOException
		{
		final SAMSequenceDictionary dict=SAMSequenceDictionaryExtractor.extractDictionary(Paths.get(inputFile));
		if(dict==null || dict.isEmpty()) return null;
		final Path dictF = support.createTmpPath(".dict");
		final BufferedWriter bw = Files.newBufferedWriter(dictF);
		final SAMSequenceDictionaryCodec codec= new SAMSequenceDictionaryCodec(bw);
		codec.encode(new SAMSequenceDictionary(
				dict.getSequences().stream().
				map(S->new SAMSequenceRecord(
						S.getSequenceName().startsWith("chr")?S.getSequenceName().substring(3):"chr"+S.getSequenceName(),
						S.getSequenceLength())).
				collect(Collectors.toList())));
		bw.flush();
		bw.close();
		return dictF;
		}
	
	@Test(dataProvider="all-vcf-files")
	public void test01(final String inputFile) 
		throws IOException
		{
		try {
		final Path dictF = prepare(inputFile);
		final Path output = support.createTmpPath(".vcf");
        Assert.assertEquals(new VcfSetSequenceDictionary().instanceMain(new String[]{
        		"-o",output.toString(),
        		"-R",dictF.toString(),
        		inputFile
        	}),0);
        support.assertIsVcf(output);
		} finally {
			support.removeTmpFiles();
		}
		}
	
	@Test(dataProvider="all-vcf-files")
	public void testHeaderOnly(final String inputFile) 
		throws IOException
		{
		try {
		final Path dictF = prepare(inputFile);
		final Path output = support.createTmpPath(".vcf");
        Assert.assertEquals(new VcfSetSequenceDictionary().instanceMain(new String[]{
        		"-o",output.toString(),
        		"-R",dictF.toString(),
        		"-ho",
        		inputFile
        	}),0);
        support.assertIsVcf(output);
		} finally {
			support.removeTmpFiles();
		}
		}

}

package com.github.lindenb.jvarkit.tools.misc;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.stream.Collectors;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceDictionaryCodec;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;

public class VcfSetSequenceDictionaryTest extends TestUtils {
	
	private File prepare(final String inputFile)  throws IOException
		{
		final SAMSequenceDictionary dict=SAMSequenceDictionaryExtractor.extractDictionary(new File(inputFile));
		if(dict==null || dict.isEmpty()) return null;
		final File dictF = super.createTmpFile(".dict");
		final BufferedWriter bw = Files.newBufferedWriter(dictF.toPath());
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
		final File dictF = prepare(inputFile);
		final File output = super.createTmpFile(".vcf");
        Assert.assertEquals(new VcfSetSequenceDictionary().instanceMain(new String[]{
        		"-o",output.getPath(),
        		"-R",dictF.getPath(),
        		inputFile
        	}),0);
        assertIsVcf(output);
		}
	
	@Test(dataProvider="all-vcf-files")
	public void testHeaderOnly(final String inputFile) 
		throws IOException
		{
		final File dictF = prepare(inputFile);
		final File output = super.createTmpFile(".vcf");
        Assert.assertEquals(new VcfSetSequenceDictionary().instanceMain(new String[]{
        		"-o",output.getPath(),
        		"-R",dictF.getPath(),
        		"-ho",
        		inputFile
        	}),0);
        assertIsVcf(output);
		}

}

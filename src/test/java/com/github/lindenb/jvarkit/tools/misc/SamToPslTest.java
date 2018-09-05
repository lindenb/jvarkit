package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.IOException;
import java.nio.file.Paths;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;

public class SamToPslTest extends TestUtils{
	@Test(dataProvider="all-sam-or-bam-files")
	public void test1(final String samFile) throws IOException {
		final SAMSequenceDictionary dict = SAMSequenceDictionaryExtractor.extractDictionary(Paths.get(samFile));
		if(dict==null || dict.isEmpty()) return;
		if(dict.getSequences().stream().
				mapToInt(L->L.getSequenceLength()).
				max().orElse(0) > 1_000_000) return;
 
		final File out = createTmpFile(".txt");
		Assert.assertEquals(
			new SamToPsl().instanceMain(newCmd().add(
					"-o",out,
					samFile
			).make()),0);
		super.assertIsNotEmpty(out);
	}
}

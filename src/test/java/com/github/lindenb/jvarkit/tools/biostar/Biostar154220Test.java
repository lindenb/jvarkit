package com.github.lindenb.jvarkit.tools.biostar;

import java.io.File;
import java.io.IOException;
import java.nio.file.Paths;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.misc.SortSamRefName;
import com.github.lindenb.jvarkit.tools.tests.TestUtils;

import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;

public class Biostar154220Test extends TestUtils {
	@Test(dataProvider="all-sam-or-bam-files")
	public void test1(final String samFile) throws IOException {
		// limit to small ref
		if(SAMSequenceDictionaryExtractor.extractDictionary(Paths.get(samFile)).
				getSequences().stream().
				mapToInt(L->L.getSequenceLength()).
				min().orElse(0) > 1_000_000) return;
		
		final File out1 = createTmpFile(".bam");
		Assert.assertEquals(
				new  SortSamRefName().instanceMain(newCmd().add(
						"-o",out1,
						samFile
				).make()),0);
		
		
		super.assertIsValidBam(out1);
		
		final File out2 = createTmpFile(".bam");
		Assert.assertEquals(
				new  Biostar154220().instanceMain(newCmd().add(
						"-o",out2,
						"-n",5,
						out1.getPath()
				).make()),0);
		
		super.assertIsValidBam(out2);
	}
}

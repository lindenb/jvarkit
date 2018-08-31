package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.IOException;
import java.nio.file.Paths;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;

public class SamAddPITest extends TestUtils {
	@Test(dataProvider="all-sam-or-bam-files")
	public void test1(final String samFile) throws IOException {
		// limit to small ref
		if(SAMSequenceDictionaryExtractor.extractDictionary(Paths.get(samFile)).
				getSequences().stream().
				mapToInt(L->L.getSequenceLength()).
				min().orElse(0) > 1_000_000) return;

		
		final File out = createTmpFile(".bam");
		final File in = addClippingToBam(new File(samFile));
		Assert.assertEquals(
			new SamAddPI().instanceMain(newCmd().add(
					"-o",out,"-w","-N","1000",
					in
			).make()),0);
		super.assertIsValidBam(out);
	}
}

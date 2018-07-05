package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.IOException;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

public class SamAddPITest extends TestUtils {
	@Test(dataProvider="all-sam-or-bam-files")
	public void test1(final String samFile) throws IOException {
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

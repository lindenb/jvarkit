package com.github.lindenb.jvarkit.tools.biostar;

import java.io.File;
import java.io.IOException;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.misc.SortSamRefName;
import com.github.lindenb.jvarkit.tools.tests.TestUtils;

public class Biostar154220Test extends TestUtils {
	@Test(dataProvider="all-sam-or-bam-files")
	public void test1(final String samFile) throws IOException {
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

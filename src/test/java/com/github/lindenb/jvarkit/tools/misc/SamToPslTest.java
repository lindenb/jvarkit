package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.IOException;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

public class SamToPslTest extends TestUtils{
	@Test(dataProvider="all-sam-or-bam-files")
	public void test1(final String samFile) throws IOException {
		final File out = createTmpFile(".txt");
		Assert.assertEquals(
			new SamToPsl().instanceMain(newCmd().add(
					"-o",out,
					samFile
			).make()),0);
		super.assertIsNotEmpty(out);
	}
}

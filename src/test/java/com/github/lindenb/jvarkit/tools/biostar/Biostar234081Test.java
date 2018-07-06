package com.github.lindenb.jvarkit.tools.biostar;

import java.io.File;
import java.io.IOException;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

public class Biostar234081Test extends TestUtils {
	@Test(dataProvider="all-sam-or-bam-files")
	public void test1(final String samFile) throws IOException {
		
		final File out = createTmpFile(".bam");
		Assert.assertEquals(
				new  Biostar234081().instanceMain(newCmd().add(
						"-o",out.getPath(),
						samFile
				).make()),0);
		
		super.assertIsValidBam(out);
	}
}

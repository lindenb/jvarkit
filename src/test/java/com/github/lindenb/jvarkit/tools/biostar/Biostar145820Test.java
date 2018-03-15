package com.github.lindenb.jvarkit.tools.biostar;

import java.io.File;
import java.io.IOException;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;


public class Biostar145820Test  extends TestUtils{
	
	@Test(dataProvider="all-sam-or-bam-files")
	public void testShuffle(final String bam) throws IOException {
		final File out = createTmpFile(".bam");
		Assert.assertEquals(
			new Biostar145820().instanceMain(newCmd().
			add("-o").add(out).
			add(bam).
			make()
			),0);
		assertIsValidBam(out);
		Assert.assertEquals(wc(out), wc(new File(bam)));
		}
	
	@Test(dataProvider="all-sam-or-bam-files")
	public void testDownSample(final String bam) throws IOException {
		final File out = createTmpFile(".bam");
		Assert.assertEquals(
			new Biostar145820().instanceMain(newCmd().
			add("-o").add(out).
			add("-n").add(10).
			add(bam).
			make()
			),0);
		assertIsValidBam(out);
		Assert.assertTrue(wc(out)<=10);
		}
	}

package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.IOException;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

public class VCFShuffleTest  extends TestUtils {
	@Test(dataProvider="all-vcf-files")
	public void test01(final String inputFile) 
		throws IOException
		{
		final File out = super.createTmpFile(".vcf");
		Assert.assertEquals(0,new VCFShuffle().instanceMain(new String[] {
			"-o",out.getPath(),
			inputFile
			}));
		assertIsVcf(out);
		Assert.assertEquals(
			super.variantStream(new File(inputFile)).count(),
			super.variantStream(out).count()
			);
		}
	}

package com.github.lindenb.jvarkit.tools.biostar;

import java.io.File;
import java.nio.file.Paths;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

public class Biostar352930Test extends TestUtils {
	@Test(dataProvider="all-sam-or-bam-files")
	public void test01(final String bampath) throws Exception
		{
		final File out = super.createTmpFile(".bam");
		final File sortedBam1= sortBamOnQueryName(Paths.get(bampath),null);
		
		Assert.assertEquals(new Biostar352930().instanceMain(new String[] {
			"-o",out.getPath(),
			sortedBam1.getPath()
			}),0);
		super.assertIsValidBam(out);
		}
	}

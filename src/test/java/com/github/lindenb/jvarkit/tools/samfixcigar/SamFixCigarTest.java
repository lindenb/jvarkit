package com.github.lindenb.jvarkit.tools.samfixcigar;

import java.io.File;
import java.io.IOException;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

public class SamFixCigarTest  extends TestUtils {
	@Test(dataProvider="all-one-bam-and-ref")
	public void test01(final String inBam,String inFasta) 
		throws IOException
		{
		final File out = super.createTmpFile(".bam");
		final SamFixCigar cmd =new SamFixCigar();
		Assert.assertEquals(cmd.instanceMain(new String[] {
			"-R",inFasta,
			"-o",out.getPath(),
			inBam
			}),0);
		assertIsValidBam(out);
		}
}

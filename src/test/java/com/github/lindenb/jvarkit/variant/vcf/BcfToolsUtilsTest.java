package com.github.lindenb.jvarkit.variant.vcf;

import java.io.IOException;
import java.nio.file.Paths;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;

public class BcfToolsUtilsTest {
	private final TestSupport support  = new TestSupport();
	@Test
	public void test1() throws IOException {
		try {
			Assert.assertTrue(BcfToolsUtils.isBcfToolsRequired(Paths.get(support.resource("toy.bcf"))));
			Assert.assertFalse(BcfToolsUtils.isBcfToolsRequired(Paths.get(support.resource("toy.vcf.gz"))));
			}
		finally
			{
			support.removeTmpFiles();
			}
		}

}

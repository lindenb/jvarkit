package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.IOException;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

public class VariantsInWindowTest  extends TestUtils {


 @Test(dataProvider="all-vcf-files")
 public void test01(final String input) throws IOException {
		final File out = super.createTmpFile(".vcf");
		Assert.assertEquals(0,new VariantsInWindow().instanceMain(new String[] {
			"-o",out.getPath(),
			input
			}));
		assertIsVcf(out);
		}
}

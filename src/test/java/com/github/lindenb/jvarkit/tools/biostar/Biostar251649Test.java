package com.github.lindenb.jvarkit.tools.biostar;

import java.io.File;
import java.io.IOException;

import org.testng.Assert;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

public class Biostar251649Test extends TestUtils {
	
	public void test01() throws IOException {
		final File out = super.createTmpFile(".vcf");
		Assert.assertEquals(
			new Biostar251649().instanceMain(newCmd().
			add("-o").add(out).
			add("-R",SRC_TEST_RESOURCE+"/rotavirus_rf.fa").
			add(SRC_TEST_RESOURCE+"/S1.vcf.gz").
			make()
			),0);
		assertIsVcf(out);
		}
}

package com.github.lindenb.jvarkit.tools.upstreamorf;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

import java.io.File;
import java.io.IOException;

import org.testng.Assert;
import org.testng.annotations.Test;

public class VcfScanUpstreamOrfTest  extends TestUtils {
	@Test
	public void test01(final String vcf) throws IOException {
		final File out = createTmpFile(".vcf");
		
		Assert.assertEquals(new VcfScanUpstreamOrf().instanceMain(new String[] {
			"-o",out.getPath(),
			"-R",SRC_TEST_RESOURCE+"/rotavirus_rf.fa",
			"-k",SRC_TEST_RESOURCE+"/rotavirus_rf.knowngenes.tsv.gz",
			SRC_TEST_RESOURCE+"/rotavirus_rf.vcf.gz"
			}),0);
		assertIsVcf(out);
		}
	}

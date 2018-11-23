package com.github.lindenb.jvarkit.tools.structvar;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

public class SamScanSplitReadsTest  extends TestUtils {
	@Test
	public void test01() throws java.io.IOException {
		final File out= super.createTmpFile(".vcf");
		Assert.assertEquals(new SamScanSplitReads().instanceMain(new String[] {
				"-o",out.getPath(),
				SRC_TEST_RESOURCE+"/S1.bam",
				SRC_TEST_RESOURCE+"/S2.bam",
				SRC_TEST_RESOURCE+"/S3.bam"
				}),0);
		super.assertIsVcf(out);
	}

}

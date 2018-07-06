package com.github.lindenb.jvarkit.tools.vcfcmp;


import java.io.File;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

public class VcfCompareCallersTest extends TestUtils {

@DataProvider(name="src01")
public Object[][] getData01() {
	return new Object[][] {
		{SRC_TEST_RESOURCE+"/rotavirus_rf.freebayes.vcf.gz",SRC_TEST_RESOURCE+"/rotavirus_rf.vcf.gz"},
		{SRC_TEST_RESOURCE+"/rotavirus_rf.vcf.gz",SRC_TEST_RESOURCE+"/rotavirus_rf.unifiedgenotyper.vcf.gz"},
		{SRC_TEST_RESOURCE+"/rotavirus_rf.freebayes.vcf.gz",SRC_TEST_RESOURCE+"/rotavirus_rf.unifiedgenotyper.vcf.gz"}
	};
	}

@Test(dataProvider="src01")
public void compareCallers(final String vcf1,final String vcf2) throws Exception {
	File outFile = super.createTmpFile(".zip");
	Assert.assertEquals(new  VcfCompareCallers().instanceMain(newCmd().add(
			"-o",outFile.getPath(),
			vcf1,vcf2
			).make()
		),0);
	assertZip(outFile);
	}
}

package com.github.lindenb.jvarkit.tools.vcfcmp;


import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;

@AlsoTest(LauncherTest.class)
public class VcfCompareCallersTest {
	private final TestSupport support =new TestSupport();

@DataProvider(name="src01")
public Object[][] getData01() {
	return new Object[][] {
		{support.resource("rotavirus_rf.freebayes.vcf.gz"),support.resource("rotavirus_rf.vcf.gz")},
		{support.resource("rotavirus_rf.vcf.gz"),support.resource("rotavirus_rf.unifiedgenotyper.vcf.gz")},
		{support.resource("rotavirus_rf.freebayes.vcf.gz"),support.resource("rotavirus_rf.unifiedgenotyper.vcf.gz")}
		};
	}

@Test(dataProvider="src01")
public void compareCallers(final String vcf1,final String vcf2) throws Exception {
	try {
	Path outFile = support.createTmpPath(".zip");
	Assert.assertEquals(new  VcfCompareCallers().instanceMain(new String[] {
			"-o",outFile.toString(),
			vcf1,vcf2
		}),0);
	support.assertZip(outFile);
	} finally {
		support.removeTmpFiles();
	}
	}
}

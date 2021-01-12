package com.github.lindenb.jvarkit.tools.gnomad;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;


@AlsoTest(LauncherTest.class)
public class VcfGnomadTest {
	private final TestSupport support = new TestSupport();

	
@DataProvider(name="src01")
public Object[][] testData01() {
	return new Object[][] {
		{support.resource("test_vcf01.vcf"),support.resource("gnomad.exomes.r2.0.1.sites.vcf.gz")},
		{support.resource("test_vcf01.vcf"),support.resource("gnomad.genomes.r2.0.1.sites.vcf.gz")}
	};
}

@Test(dataProvider="src01")
public void test01(final String vcfpath,final String gnomad) throws IOException {
	try {
		final Path vcfOut = support.createTmpPath(".vcf");
		
		Assert.assertEquals(new VcfGnomad().instanceMain(new String[]{
				"-o",vcfOut.toString(),
				"-g",gnomad,
				vcfpath
				}),0);
		support.assertIsVcf(vcfOut);
		}
	finally
		{
		support.removeTmpFiles();
		}
	}



}

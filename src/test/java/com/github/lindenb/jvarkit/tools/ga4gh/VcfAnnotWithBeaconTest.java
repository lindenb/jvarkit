package com.github.lindenb.jvarkit.tools.ga4gh;

import java.io.File;
import java.io.IOException;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

public class VcfAnnotWithBeaconTest extends TestUtils {

	@DataProvider(name = "src1")
	public Object[][] createData1() {
		return new Object[][]{
			{SRC_TEST_RESOURCE+"/test_vcf01.vcf"},
			{SRC_TEST_RESOURCE+"/gnomad.genomes.r2.0.1.sites.1.vcf.gz"},
			{SRC_TEST_RESOURCE+"/gnomad.exomes.r2.0.1.sites.1.vcf.gz"},
			{SRC_TEST_RESOURCE+"/ExAC.r1.sites.vep.vcf.gz"}
			};
		}
	
	@Test(dataProvider="src1")
	public void test1(final String vcf) throws IOException {
		final File out = createTmpFile(".vcf");
		Assert.assertEquals(
			new VcfAnnotWithBeacon().instanceMain(newCmd().
			add("-o").add(out).
			add("--cert").
			add(vcf).
			make()
			),0);
		assertIsVcf(out);
		}
}

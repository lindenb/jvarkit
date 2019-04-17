package com.github.lindenb.jvarkit.tools.ga4gh;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;

public class VcfAnnotWithBeaconTest  {
	private final TestSupport support = new TestSupport();

	@DataProvider(name = "src1")
	public Object[][] createData1() {
		return new Object[][]{
			{support.resource("test_vcf01.vcf")},
			{support.resource("gnomad.genomes.r2.0.1.sites.1.vcf.gz")},
			{support.resource("gnomad.exomes.r2.0.1.sites.vcf.gz")},
			{support.resource("/ExAC.r1.sites.vep.vcf.gz")}
			};
		}
	
	@Test(dataProvider="src1",enabled=false)
	public void test1(final String vcf) throws IOException {
			try {
			final Path out = support.createTmpPath(".vcf");
			Assert.assertEquals(
				new VcfAnnotWithBeacon().instanceMain(new String[] {
				"-o",out.toString(),
				"--cert",
				"--tee",
				vcf}
				),0);
			support.assertIsVcf(out);
			}
		finally {
			support.removeTmpFiles();
			}
		}
}

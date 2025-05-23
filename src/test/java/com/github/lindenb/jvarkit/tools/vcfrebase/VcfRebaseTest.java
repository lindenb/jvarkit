package com.github.lindenb.jvarkit.tools.vcfrebase;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;


import com.github.lindenb.jvarkit.tools.tests.TestSupport;


public class VcfRebaseTest {

	@DataProvider(name = "src1")
	public Object[][] createData1() {
		return new Object[][]{
			{"S1.vcf.gz","rotavirus_rf.fa"},
			{"toy.vcf.gz","toy.fa"}
			};
		}
	@Test(dataProvider="src1")
	public void test1(final String vcf,final String ref) throws IOException {
		 final TestSupport support = new TestSupport();
		try {
			final Path out = support.createTmpPath(".vcf");
			Assert.assertEquals(
				new VcfRebase().instanceMain(new String[] {
				"-o",out.toString(),
				"-R",support.resource(ref),
				"-E","EcoRI",
				"-E","BamHI",
				support.resource(vcf)}
				),0);
			support.assertIsVcf(out);
			}
		finally {
			support.removeTmpFiles();
			}
		}
	@Test
	public void test2() throws IOException {
		final TestSupport support = new TestSupport();
		try {
			final Path ref=	support.getGRCh37Path().orElse(null);
			if(ref==null) return;
			final Path out = support.createTmpPath(".vcf");
			Assert.assertEquals(
				new VcfRebase().instanceMain(new String[] {
				"-o",out.toString(),
				"-R",ref.toString(),
				"-E","EcoRI",
				"-E","BamHI",
				support.resource("gnomad.exomes.r2.0.1.sites.vcf.gz")
				}),0);
			support.assertIsVcf(out);
			}
		finally {
			support.removeTmpFiles();
			}
		}
}

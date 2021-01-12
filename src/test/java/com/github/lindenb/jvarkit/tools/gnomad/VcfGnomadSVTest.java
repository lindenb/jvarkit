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
public class VcfGnomadSVTest {
	private final TestSupport support = new TestSupport();

	
@DataProvider(name="src01")
public Object[] testData01() {
	return new Object[] {
			support.resource("manta.B00GWGD.vcf.gz"),
			support.resource("manta.B00GWIU.vcf.gz"),
			support.resource("manta.B00I9CJ.vcf.gz"),
			support.resource("manta.D000Q1R.vcf.gz")	
		};
	}

@Test(dataProvider="src01")
public void test01(final String vcfpath) throws IOException {
	try {
		final Path vcfOut = support.createTmpPath(".vcf");
		
		Assert.assertEquals(new VcfGnomadSV().instanceMain(new String[]{
				"-o",vcfOut.toString(),
				"-g",support.resource("gnomad_v2_sv.sites.vcf.gz"),
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

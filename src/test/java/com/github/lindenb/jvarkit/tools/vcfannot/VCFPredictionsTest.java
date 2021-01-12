package com.github.lindenb.jvarkit.tools.vcfannot;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.bio.structure.GtfReaderTest;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;

@AlsoTest({LauncherTest.class,GtfReaderTest.class})
public class VCFPredictionsTest {
	private final TestSupport support =new TestSupport();
	@DataProvider(name="data1")
	public Object[] data01() {
		return new Object[] {
				support.resource("test_vcf01.vcf"),
				support.resource("ExAC.r1.sites.vep.vcf.gz"),
				support.resource("gnomad.exomes.r2.0.1.sites.vcf.gz"),
				support.resource("gnomad.genomes.r2.0.1.sites.1.vcf.gz"),
				support.resource("manta.B00GWGD.vcf.gz")
				
		};
	}
	
	@Test(dataProvider="data1")
	public void test01(String vcfin) throws IOException
		{
		try {
			final Path ref = support.getGRCh37Path().orElse(null);
			if(ref==null) return;
			Path out = support.createTmpPath(".vcf");
			Assert.assertEquals(new VCFPredictions().instanceMain(new String[] {
				"-o",out.toString(),
				"--gtf",support.resource("Homo_sapiens.GRCh37.87.gtf.gz"),
				"-R",ref.toString(),
				vcfin
				}),0);
			support.assertIsVcf(out);
			}
		finally
			{
			support.removeTmpFiles();
			}
		}
}

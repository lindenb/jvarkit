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

public class SVPredictionsTest {
	private final TestSupport support =new TestSupport();

	@DataProvider(name="data1")
	public Object[] data01() {
		return new Object[] {
				support.resource("manta.B00GWGD.vcf.gz"),
				support.resource("manta.B00GWIU.vcf.gz"),
				support.resource("manta.B00I9CJ.vcf.gz"),
				support.resource("manta.D000Q1R.vcf.gz")				
				
		};
	}
	
	@Test(dataProvider="data1")
	public void test01(String vcfin) throws IOException
		{
		try {
			Path out = support.createTmpPath(".vcf");
			Assert.assertEquals(new SVPredictions().instanceMain(new String[] {
				"-o",out.toString(),
				"--gtf",support.resource("Homo_sapiens.GRCh37.87.gtf.gz"),
				vcfin
				}),0);
			support.assertIsVcf(out);
			}
		finally
			{
			support.removeTmpFiles();
			}
		}
	@Test(dataProvider="data1")
	public void test02(String vcfin) throws IOException
		{
		try {
			Path out = support.createTmpPath(".vcf");
			Assert.assertEquals(new SVPredictions().instanceMain(new String[] {
				"-o",out.toString(),
				"--where","intergenic,exon,intron",
				"--filter","XXXX",
				"--remove-attribute",
				"--gtf",support.resource("Homo_sapiens.GRCh37.87.gtf.gz"),
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

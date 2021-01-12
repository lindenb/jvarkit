package com.github.lindenb.jvarkit.tools.vcfannot;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.bio.GranthamScoreTest;
import com.github.lindenb.jvarkit.util.bio.structure.GtfReaderTest;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;

@AlsoTest({LauncherTest.class,GtfReaderTest.class,GranthamScoreTest.class})
public class VCFCombineTwoSnvsTest {
	private final TestSupport support =new TestSupport();
	
	@Test
	public void test01() throws IOException
		{
		try {
			Path out = support.createTmpPath(".vcf");
			Assert.assertEquals(new VCFCombineTwoSnvs().instanceMain(new String[] {
				"-R",support.resource("rotavirus_rf.fa"),
				"-o",out.toString(),
				"--gtf",support.resource("Homo_sapiens.GRCh37.87.gtf.gz"),
				"-B",support.resource("S1.bam"),
				support.resource("S1.vcf.gz")
				}),0);
			support.assertIsVcf(out);
			}
		finally
			{
			support.removeTmpFiles();
			}
		}
	
}

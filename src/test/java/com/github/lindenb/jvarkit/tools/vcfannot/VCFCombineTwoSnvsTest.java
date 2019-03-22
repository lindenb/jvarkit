package com.github.lindenb.jvarkit.tools.vcfannot;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.bio.GranthamScoreTest;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;
import com.github.lindenb.jvarkit.util.ucsc.KnownGeneTest;

@AlsoTest({LauncherTest.class,KnownGeneTest.class,GranthamScoreTest.class})
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
				"-k",support.resource("rotavirus_rf.knowngenes.tsv.gz"),
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

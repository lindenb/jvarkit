package com.github.lindenb.jvarkit.tools.upstreamorf;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.io.ArchiveFactoryTest;
import com.github.lindenb.jvarkit.tools.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.AlgorithmsTest;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;
import com.github.lindenb.jvarkit.util.ucsc.KnownGeneTest;
import com.github.lindenb.jvarkit.util.vcf.VCFUtilsTest;

@AlsoTest(value= {LauncherTest.class,AlgorithmsTest.class,ArchiveFactoryTest.class,KnownGeneTest.class,VCFUtilsTest.class})
public class VcfScanUpstreamOrfTest {
	private final TestSupport support =new TestSupport();
	@Test
	public void testBasicVcf() throws IOException {
		try {
			Path out = support.createTmpPath(".vcf");
			
			Assert.assertEquals(new VcfScanUpstreamOrf().instanceMain(new String[] {
				"-o",out.toString(),
				"-R",support.resource("rotavirus_rf.fa"),
				"-k",support.resource("rotavirus_rf.knowngenes.tsv.gz"),
				support.resource("rotavirus_rf.vcf.gz")
				}),0);
			support.assertIsVcf(out);
			}
		finally
			{
			support.removeTmpFiles();
			}
		}
	}

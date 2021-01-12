package com.github.lindenb.jvarkit.tools.upstreamorf;

import java.io.IOException;
import java.nio.file.Path;
import java.util.Optional;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.io.ArchiveFactoryTest;
import com.github.lindenb.jvarkit.tests.AlsoTest;
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
			Optional<Path> ref=support.getGRCh37Path();
			if(!ref.isPresent()) return;
			
			Path out = support.createTmpPath(".vcf");
			Assert.assertEquals(new VcfScanUpstreamOrf().instanceMain(new String[] {
				"-o",out.toString(),
				"-R",ref.get().toString(),
				"--gtf",support.resource("Homo_sapiens.GRCh37.87.gtf.gz"),
				support.resource("test_vcf01.vcf")
				}),0);
			support.assertIsVcf(out);
			}
		finally
			{
			support.removeTmpFiles();
			}
		}
	}

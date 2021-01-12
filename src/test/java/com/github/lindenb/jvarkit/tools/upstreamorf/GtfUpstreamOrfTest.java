package com.github.lindenb.jvarkit.tools.upstreamorf;

import java.io.IOException;
import java.nio.file.Path;
import java.util.Optional;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.bio.KozakSequenceTest;
import com.github.lindenb.jvarkit.util.bio.structure.GtfReader;
import com.github.lindenb.jvarkit.util.bio.structure.GtfReaderTest;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;
import com.github.lindenb.jvarkit.util.ucsc.KnownGeneTest;
import com.github.lindenb.jvarkit.util.vcf.VCFUtilsTest;


@AlsoTest(value= {LauncherTest.class,KozakSequenceTest.class,GtfReaderTest.class,KnownGeneTest.class,VCFUtilsTest.class})
public class GtfUpstreamOrfTest {
	private final TestSupport support =new TestSupport();
	@Test
	public void testBasicVcf() throws IOException {
		try {
			Optional<Path> ref = support.getGRCh37Path();
			if(!ref.isPresent()) return;
			
			Path out = support.createTmpPath(".gtf");
			
			Assert.assertEquals(new GtfUpstreamOrf().instanceMain(new String[] {
				"-o",out.toString(),
				"-R",ref.get().toString(),
				support.resource("Homo_sapiens.GRCh37.87.gtf.gz")
				}),0);
			
			GtfReader gf = new GtfReader(out);
			gf.getAllGenes();
			gf.close();
			
			}
		finally
			{
			support.removeTmpFiles();
			}
		}
	}

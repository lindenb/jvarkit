package com.github.lindenb.jvarkit.tools.liftover;

import java.io.IOException;
import java.net.URL;
import java.nio.file.Path;
import java.util.zip.GZIPInputStream;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.tools.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;


@AlsoTest(LauncherTest.class)
public class VcfLiftOverTest {
	final TestSupport support= new TestSupport();
	
	@Test
	public void test01() throws IOException {
		try {
			final String url="http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz";
			Path chain = support.createTmpPath(".chain");
			GZIPInputStream in=new GZIPInputStream(new URL(url).openStream());
			IOUtils.copyTo(in,chain);
			in.close();
			support.assertIsNotEmpty(chain);
			final Path out = support.createTmpPath(".vcf");

			
			Assert.assertEquals(new VcfLiftOver().instanceMain(new String[] {
					"-o",out.toString(),
					"--chain",chain.toString(),
					"-R",support.resource("rotavirus_rf.fa"),
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

package com.github.lindenb.jvarkit.tools.misc;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;
import com.github.lindenb.jvarkit.variant.variantcontext.AttributeCleanerTest;


@AlsoTest({LauncherTest.class,AttributeCleanerTest.class})
public class VcfGtfSplitterTest {
	private final TestSupport support = new TestSupport();

	@Test
	public void test01() throws IOException {
		try {
			Path tmp = support.createTmpPath(".zip");
			Path manifest = support.createTmpPath(".mf");
			
			Assert.assertEquals(new VcfGtfSplitter().instanceMain(new String[] {
					"-m",manifest.toString(),
					"--gtf",support.resource("Homo_sapiens.GRCh37.87.gtf.gz"),
					"--index",
					"-o",tmp.toString(),
					"--xannotate","INFO/ANN",
					support.resource("test_vcf01.vcf")
					}),0);
			
			support.assertZip(tmp);
			support.assertIsBed(manifest);
			support.assertTsvTableIsConsitent(manifest, null);
			}
		finally {
			support.removeTmpFiles();
		}
	}
	
}

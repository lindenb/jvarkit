package com.github.lindenb.jvarkit.tools.misc;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;

public class VcfGeneSplitterTest {
	private final TestSupport support = new TestSupport();

	public void test01() throws IOException {
		try {
			Path tmp = support.createTmpPath(".zip");
			Path manifest = support.createTmpPath(".mf");
			
			Assert.assertEquals(new VcfGeneSplitter().instanceMain(new String[] {
					"-m",manifest.toString(),
					"-o",tmp.toString(),
					support.resource("rotavirus_rf.ann.vcf.gz")
					}),0);
			
			support.assertZip(tmp);
			support.assertIsBed(manifest);
			}
		finally {
			support.removeTmpFiles();
		}
	}
}

package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;

import htsjdk.samtools.util.IOUtil;

@AlsoTest(LauncherTest.class)
public class VcfGeneSplitterTest {
	private final TestSupport support = new TestSupport();

	@Test
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
			support.assertTsvTableIsConsitent(manifest, null);
			}
		finally {
			support.removeTmpFiles();
		}
	}
	
	@Test
	public void test02() throws IOException {
		File tmp=null;
		try {
			tmp = IOUtil.createTempDir("tmp.", ".dir");
			Path manifest = support.createTmpPath(".mf");
			
			Assert.assertEquals(new VcfGeneSplitter().instanceMain(new String[] {
					"-m",manifest.toString(),
					"-o",tmp.toString(),
					support.resource("rotavirus_rf.ann.vcf.gz")
					}),0);
			
			support.assertIsBed(manifest);
			support.assertTsvTableIsConsitent(manifest, null);
			}
		finally {
			if(tmp!=null) IOUtil.deleteDirectoryTree(tmp);
			support.removeTmpFiles();
		}
	}
}

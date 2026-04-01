package com.github.lindenb.jvarkit.tools.vcfsplitgene;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.Test;


import com.github.lindenb.jvarkit.tools.tests.TestSupport;

import htsjdk.samtools.util.IOUtil;


public class VcfGeneSplitterTest {

	@Test
	public void inZip() throws IOException {
		final TestSupport support = new TestSupport();
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
	public void inDirectory() throws IOException {
		final TestSupport support = new TestSupport();
		Path tmp=null;
		try {
			tmp = IOUtil.createTempDir("tmp.");
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
			if(tmp!=null) IOUtil.deleteDirectoryTree(tmp.toFile());
			support.removeTmpFiles();
		}
	}
}

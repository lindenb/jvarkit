package com.github.lindenb.jvarkit.tools.vcfsplit;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.Test;


import com.github.lindenb.jvarkit.tools.tests.TestSupport;

import htsjdk.samtools.util.IOUtil;


public class VcfSplitNVariantsTest {
	
	private void vcfSplit(String key1,String value1) throws IOException {
		 final TestSupport support = new TestSupport();
		 Path tmpDir  = null;
		try {
			tmpDir = IOUtil.createTempDir("tmp");
			Path manifest = support.createTmpPath(".mf");
			
			
			Assert.assertEquals(new VcfSplitNVariants().instanceMain(new String[] {
					"-m",manifest.toString(),
					"--force",
					key1,value1,
					"-o",tmpDir.toString()+File.separatorChar+"TMP",
					support.resource("rotavirus_rf.vcf.gz")
					}),0);			
			}
		finally {
			support.removeTmpFiles();
			if(tmpDir!=null) IOUtil.deleteDirectoryTree(tmpDir.toFile());
		}
	}
	
	@Test
	public void testDistance() throws IOException {
		vcfSplit("--distance","1000");
	}
	
	@Test
	public void testVcfCount() throws IOException {
		vcfSplit("--vcf-count","3");
	}
	@Test
	public void testVariantsCount() throws IOException {
		vcfSplit("--vc-count","100");
	}
}

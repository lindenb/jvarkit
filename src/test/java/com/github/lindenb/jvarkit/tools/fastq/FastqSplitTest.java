package com.github.lindenb.jvarkit.tools.fastq;

import java.io.File;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;

import htsjdk.samtools.util.IOUtil;

@AlsoTest(LauncherTest.class)
public class FastqSplitTest {
	private final TestSupport support = new TestSupport();

	@Test
	public void test01() {
		try {
			File tmp = IOUtils.createTempDir("split", ".tmp", IOUtils.getDefaultTmpDir());
			Assert.assertEquals(new FastqSplit().instanceMain(new String[] {
					"-o",tmp.toString()+File.separator+"jeter.split.__TOKEN__.fastq.gz",
					"-n","100",
					"--validate",
					support.resource("S1.R1.fq.gz"),
					support.resource("S1.R2.fq.gz"),
					}),0);
			final File files[]=tmp.listFiles();
			Assert.assertTrue(files.length>0);
			for(File f:files) {
				support.assertIsFastq(f.toPath());
				}
			IOUtil.deleteDirectoryTree(tmp);
		} finally {
			support.removeTmpFiles();
		}
	}
	
	
	@Test
	public void test02() {
		try {
			File tmp = IOUtils.createTempDir("split", ".tmp", IOUtils.getDefaultTmpDir());
			Assert.assertEquals(new FastqSplit().instanceMain(new String[] {
					"-o",tmp.toString()+File.separator+"jeter.split.__TOKEN__.fastq.gz",
					"-S","5",
					"--validate",
					support.resource("S1.R1.fq.gz"),
					support.resource("S1.R2.fq.gz"),
					}),0);
			final File files[]=tmp.listFiles();
			Assert.assertTrue(files.length==5*2);
			for(File f:files) {
				support.assertIsFastq(f.toPath());
				}
			IOUtil.deleteDirectoryTree(tmp);
		} finally {
			support.removeTmpFiles();
		}
	}

	@Test
	public void test03() {
		try {
			File tmp = IOUtils.createTempDir("split", ".tmp", IOUtils.getDefaultTmpDir());
			Assert.assertEquals(new FastqSplit().instanceMain(new String[] {
					"-o",tmp.toString()+File.separator+"jeter.split.__TOKEN__.fastq.gz",
					"-S","5",
					"-oi",
					"--validate",
					support.resource("S1.R1.fq.gz"),
					support.resource("S1.R2.fq.gz"),
					}),0);
			final File files[]=tmp.listFiles();
			Assert.assertTrue(files.length==5);
			for(File f:files) {
				support.assertIsFastq(f.toPath());
				}
			IOUtil.deleteDirectoryTree(tmp);
		} finally {
			support.removeTmpFiles();
		}
	}

	
}

package com.github.lindenb.jvarkit.tools.cmpbams;

import java.io.File;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;

@AlsoTest(LauncherTest.class)
public class CompareBamsTest{
	private final TestSupport support = new TestSupport();

		
	@Test
	public void test1() throws Exception
		{
		try {
		final Path out = support.createTmpPath(".txt");
		final String bam = support.resource("toy.bam");
		final Path tmpBam = support.createTmpPath(".bam");
		try(InputStream in=Files.newInputStream(Paths.get(bam))) {
			IOUtils.copyTo(in, tmpBam);
		}
		
		Assert.assertEquals(new CompareBams().instanceMain(new String[] {
			"-o",out.toString(),
			tmpBam.toString(),
			bam,
			}),0);
		} finally {
			support.removeTmpFiles();
		}
		}
}

package com.github.lindenb.jvarkit.tools.burden;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;

@AlsoTest(LauncherTest.class)
public class VcfBurdenSlidingWindowTest {
	private final TestSupport support = new TestSupport();

	@Test
	public void test01() throws IOException {
		try {
			final Path output = support.createTmpPath(".tsv");
			
			Assert.assertEquals(new VcfBurdenSlidingWindow().instanceMain(new String[] {
	        		"-o",output.toString(),
	        		"--pedigree",support.resource("test_vcf01.ped"),
	        		support.resource("test_vcf01.vcf")
					}),0);
			support.assertTsvTableIsConsitent(output, null);
		} finally {
			support.removeTmpFiles();
		}
	}
}

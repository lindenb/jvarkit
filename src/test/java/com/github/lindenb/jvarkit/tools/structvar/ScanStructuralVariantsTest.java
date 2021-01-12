package com.github.lindenb.jvarkit.tools.structvar;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;

@AlsoTest(LauncherTest.class)
public class ScanStructuralVariantsTest {
	final TestSupport support = new TestSupport();

	@Test
	public void test01() throws IOException
		{
		try {
		final Path ctrls= support.createTmpPath(".list");
		Files.write(ctrls, Arrays.asList(
				support.resource("manta.B00GWGD.vcf.gz"),
				support.resource("manta.B00GWIU.vcf.gz"),
				support.resource("manta.B00I9CJ.vcf.gz")
				));
		
		final Path out = support.createTmpPath(".vcf");
		Assert.assertEquals(new ScanStructuralVariants().instanceMain(new String[] {
				"-o",out.toString(),
				"--controls",ctrls.toString(),
				support.resource("manta.D000Q1R.vcf.gz")
			}),0
			);
		support.assertIsVcf(out);
		} finally {
			support.removeTmpFiles();
			}
		}
}

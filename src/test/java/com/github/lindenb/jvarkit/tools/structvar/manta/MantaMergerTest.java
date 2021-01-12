package com.github.lindenb.jvarkit.tools.structvar.manta;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;
import com.github.lindenb.jvarkit.variant.sv.StructuralVariantComparatorTest;

@AlsoTest({LauncherTest.class,StructuralVariantComparatorTest.class})
public class MantaMergerTest {
	private final TestSupport support = new TestSupport();
	@Test
	public void test01() throws IOException
		{
		try {
			final Path tmp = support.createTmpPath(".txt");
			Files.write(tmp, Arrays.asList(
					support.resource("manta.B00GWGD.vcf.gz"),
					support.resource("manta.B00GWIU.vcf.gz"),
					support.resource("manta.B00I9CJ.vcf.gz"),
					support.resource("manta.D000Q1R.vcf.gz")
					));
			
			final Path out = support.createTmpPath(".vcf");
			Assert.assertEquals(new MantaMerger().instanceMain(new String[] {
				"-o",out.toString(),
				tmp.toString()}),0);
			support.assertIsVcf(out);
			}
		catch(final Throwable err) {
			support.removeTmpFiles();
			}
		}
}
package com.github.lindenb.jvarkit.tools.structvar;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;

@AlsoTest(LauncherTest.class)
public class VcfGapFrequentTest{
	
	private final TestSupport support = new TestSupport();
	@DataProvider(name = "src1")
	public Object[][] createData1() {
		return support.toArrayArray(support.
				allVcfOrBcf().
				map(F->new Object[] {F})
				)
				;
		}
	
	
	@Test(dataProvider="src1")
	public void test01(final String input) 
		throws IOException
		{
		try {
		final Path out = support.createTmpPath(".bed");
		Assert.assertEquals(new VcfGapFrequent().instanceMain(new String[] {
			"-o",out.toString(),
			"--fields","AF_NFE",
			"-D",support.resource("gnomad.exomes.r2.0.1.sites.vcf.gz"),
			"-D",support.resource("gnomad.genomes.r2.0.1.sites.1.vcf.gz"),
			input
			}),0);
		Assert.assertTrue(Files.exists(out));
		} finally {
			support.removeTmpFiles();
		}
		}
	}

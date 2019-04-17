package com.github.lindenb.jvarkit.tools.vcfrebase;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.bio.RebaseTest;

@AlsoTest(RebaseTest.class)
public class VcfRebaseTest {
	private final TestSupport support = new TestSupport();

	@DataProvider(name = "src1")
	public Object[][] createData1() {
		return new Object[][]{
			{support.resource("S1.vcf.gz"),support.resource("rotavirus_rf.fa")},
			{support.resource("toy.vcf.gz"),support.resource("toy.fa")}
			};
		}
	@Test(dataProvider="src1")
	public void test1(final String vcf,final String ref) throws IOException {
			try {
			final Path out = support.createTmpPath(".vcf");
			Assert.assertEquals(
				new VcfRebase().instanceMain(new String[] {
				"-o",out.toString(),
				"-R",ref,
				"-E","EcoRI",
				"-E","BamHI",
				vcf}
				),0);
			support.assertIsVcf(out);
			}
		finally {
			support.removeTmpFiles();
			}
		}
	
}

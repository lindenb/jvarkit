package com.github.lindenb.jvarkit.tools.misc;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.vcf.VCFUtilsTest;

@AlsoTest(VCFUtilsTest.class)
public class VCFPolyXTest {
	
	private final TestSupport support = new TestSupport();

	
	
		@DataProvider(name = "src1")
		public Object[][] createData1() {
			return new Object[][]{
				{support.resource("toy.vcf.gz"),support.resource("toy.fa"),1},
				{support.resource("S1.vcf.gz"),support.resource("rotavirus_rf.fa"),2}
				};
		}
		
	@Test(dataProvider="src1")
	public void test1(final String vcf,final String ref,int n) throws IOException {
		try {
			final Path out = support.createTmpPath(".vcf");
			Assert.assertEquals( new VCFPolyX().instanceMain(new String[] {
						"-o",out.toString(),
						"-R",ref,
						"-n","1",
						vcf
				}),0);
			support.assertIsVcf(out);
			}
		finally {
			support.removeTmpFiles();
			}
		}
}

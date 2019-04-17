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
public class VariantsInWindowTest  {
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
 public void test01(final String input) throws IOException {
	 try {
		final Path out = support.createTmpPath(".vcf");
		Assert.assertEquals(0,new VariantsInWindow().instanceMain(new String[] {
			"-o",out.toString(),
			input
			}));
		support.assertIsVcf(out);
		}
	finally {
		 support.removeTmpFiles();
		 }
 	}
}

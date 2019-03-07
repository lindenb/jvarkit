package com.github.lindenb.jvarkit.tools.misc;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.io.IOUtilsTest;
import com.github.lindenb.jvarkit.tools.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;

@AlsoTest(IOUtilsTest.class)
public class SortSamRefNameTest  {
	private final TestSupport support = new TestSupport();
	@DataProvider(name = "src1")
	public Object[][] createData1() {
		return support.toArrayArray(support.
				allSamOrBams().
				map(F->new Object[] {F})
				)
				;
		}
	
	@Test(dataProvider="src1")
	public void test1(final String samFile) throws IOException {
		try {
		final Path out = support.createTmpPath(".bam");
		Assert.assertEquals(
				new  SortSamRefName().instanceMain(new String[] {
					"-o",out.toString(),
					samFile
				}),0);
		support.assertIsValidBam(out);
		} finally {
			support.removeTmpFiles();
		}
	}
}

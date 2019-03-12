package com.github.lindenb.jvarkit.tools.structvar;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;

@AlsoTest(LauncherTest.class)
public class SamExtractClipTest {
	private final TestSupport support = new TestSupport();

	@DataProvider(name = "src1")
	public Object[][] createData1() {
		return support.toArrayArray(
				support.allSamOrBams().
				map(F->new Object[] {F})
				);
		}
	
	@Test(dataProvider="src1")
	public void test01(final String inBam) 
		throws IOException
		{
		try {
			final Path out = support.createTmpPath(".fastq");
			Assert.assertEquals(new SamExtractClip().instanceMain(new String[] {
				"-o",out.toString(),
				inBam
				}),0);
			support.assertIsFastq(out);
			}
		finally {
			support.removeTmpFiles();
			}
		}
}

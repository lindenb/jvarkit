package com.github.lindenb.jvarkit.tools.misc;

import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;

@AlsoTest(LauncherTest.class)
public class SamClipIndelFractionTest{
	
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
		final Path out = support.createTmpPath(".txt");
		final Path in = support.addClippingToBam(Paths.get(samFile));
		Assert.assertEquals(
			new SamClipIndelFraction().instanceMain(new String[] {
					"-o",out.toString(),
					in.toString()
			}),0);
		support.assertIsNotEmpty(out);
		}
		finally {
			support.removeTmpFiles();
		}
	}
}

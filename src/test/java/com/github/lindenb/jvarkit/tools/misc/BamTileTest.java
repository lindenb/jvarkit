package com.github.lindenb.jvarkit.tools.misc;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;


import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;


public class BamTileTest {
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
		Assert.assertEquals(
			new SamClipIndelFraction().instanceMain(new String[] {
					"-o",out.toString(),
					samFile
				}),0);
		support.assertIsNotEmpty(out);
		} finally {
			support.removeTmpFiles();
		}
	}
}

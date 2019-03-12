package com.github.lindenb.jvarkit.tools.sam2tsv;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;

@AlsoTest(LauncherTest.class)
public class PrettySamTest {
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
	public void test01(final String inBam,String inFasta) 
		throws IOException
		{
		try {
			final Path out = support.createTmpPath(".txt");
			final PrettySam cmd =new PrettySam();
			Assert.assertEquals(cmd.instanceMain(new String[] {
				"-R",inFasta,
				"-o",out.toString(),
				inBam
				}),0);
			support.assertIsNotEmpty(out);
			}
		finally {
			this.support.removeTmpFiles();
			}
		}
}

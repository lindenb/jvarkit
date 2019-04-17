package com.github.lindenb.jvarkit.tools.misc;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;

@AlsoTest(LauncherTest.class)
public class VcfToTableTest  {

	private final TestSupport support = new TestSupport();

	@DataProvider(name="src01")
	public Object[][] testData01() {
			return support.toArrayArray(
					support.allVcfOrBcf().
					map(S->new Object[] {S})
					);
			}

	@Test(dataProvider="src01")
	public void test01(final String inputFile) 
		throws IOException
		{
		try {
			final Path out = support.createTmpPath(".vcf");
			final VcfToTable cmd =new VcfToTable();
			Assert.assertEquals(0,cmd.instanceMain(new String[] {
				"-o",out.toString(),
				inputFile
				}));
			}
		finally {
			support.removeTmpFiles();
			}
		}
}

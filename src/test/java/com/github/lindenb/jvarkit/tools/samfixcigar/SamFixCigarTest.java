package com.github.lindenb.jvarkit.tools.samfixcigar;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;

@AlsoTest(LauncherTest.class)
public class SamFixCigarTest {
	
	private final TestSupport support = new TestSupport();
	
	@DataProvider(name = "src1")
	public Object[][] createData1() {
		return new Object[][] {
			{support.resource("toy.bam"),support.resource("toy.fa")},
			{support.resource("S1.bam"),support.resource("rotavirus_rf.fa")}
		};
	}
	
	@Test(dataProvider="src1")
	public void test01(final String inBam,final String inFasta) 
		throws IOException
		{
		try {
			final Path out = support.createTmpPath(".bam");
			Assert.assertEquals(new SamFixCigar().instanceMain(new String[] {
				"-R",inFasta,
				"-o",out.toString(),
				inBam
				}),0);
			support.assertIsValidBam(out);
		} finally
			{
			support.removeTmpFiles();
			}
		}
}

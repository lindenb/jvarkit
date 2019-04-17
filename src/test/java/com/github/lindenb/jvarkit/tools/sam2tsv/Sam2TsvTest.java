package com.github.lindenb.jvarkit.tools.sam2tsv;

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
public class Sam2TsvTest {
	private final TestSupport support = new TestSupport();

	
	@DataProvider(name = "src1")
	public Object[][] createData1() {
		return new Object[][] {
			{support.resource("toy.bam"),support.resource("toy.fa")},
			{support.resource("S1.bam"),support.resource("rotavirus_rf.fa")}
		};
	}
	
	
	@Test(dataProvider="src1")
	public void test01(final String inBam,String inFasta) 
		throws IOException
		{
		try {
			final Path out = support.createTmpPath(".tsv");
			final Path bam2 = support.addClippingToBam(Paths.get(inBam));
			
			Assert.assertEquals(new Sam2Tsv().instanceMain(new String[] {
				"-R",inFasta,
				"-o",out.toString(),
				bam2.toString()
				}),0);
			support.assertTsvTableIsConsitent(out,null);
		} finally {
			support.removeTmpFiles();
		}
		}
}

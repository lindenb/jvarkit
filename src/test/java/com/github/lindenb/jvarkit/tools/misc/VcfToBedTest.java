package com.github.lindenb.jvarkit.tools.misc;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.samtools.util.IntervalExtender;
import com.github.lindenb.jvarkit.samtools.util.IntervalExtenderTest;
import com.github.lindenb.jvarkit.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;

@AlsoTest({LauncherTest.class,IntervalExtenderTest.class})
public class VcfToBedTest  {

	private final TestSupport support = new TestSupport();

	@DataProvider(name="src1")
	public Object[][] testData01() {
			return new Object[][] {
				{support.resource("manta.B00GWIU.vcf.gz")},
				{support.resource("toy.vcf.gz")},
				{support.resource("rotavirus_rf.ann.vcf.gz")}
				};
				
			}
	
@Test(dataProvider="src1")
public void test01(final String inputFile) 
	throws IOException
	{
	try {
		final Path out = support.createTmpPath(".bed");
		Assert.assertEquals(0,new VcfToBed().instanceMain(new String[] {
			inputFile
			}));
		support.assertIsBed(out);
		}
	finally {
		support.removeTmpFiles();
		}
	}
@Test(dataProvider="src1")
public void test02(final String inputFile) 
	throws IOException
	{
	try {
		final Path out = support.createTmpPath(".bed");
		Assert.assertEquals(0,new VcfToBed().instanceMain(new String[] {
			"-header",
			"--slop","200%",
			inputFile
			}));
		support.assertIsBed(out);
		}
	finally {
		support.removeTmpFiles();
		}
	}

}

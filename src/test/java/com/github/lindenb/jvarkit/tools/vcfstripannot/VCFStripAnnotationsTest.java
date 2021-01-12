package com.github.lindenb.jvarkit.tools.vcfstripannot;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;
import com.github.lindenb.jvarkit.variant.variantcontext.AttributeCleanerTest;

@AlsoTest({LauncherTest.class,AttributeCleanerTest.class})
public class VCFStripAnnotationsTest {
	private final TestSupport support = new TestSupport();
	
	@DataProvider(name = "src1")
	public Object[][] createData1() {
		String vcf1 = support.resource("rotavirus_rf.ann.vcf.gz");
		return new Object[][] {
			{vcf1,"INFO"},
			{vcf1,"FILTER"},
			{vcf1,"ID"},
			{vcf1,"QUAL"},
			{vcf1,"FORMAT"},
			{vcf1,"INFO/AC"},
			{vcf1,"^INFO/AC,INFO/AF"}
		};
		}
	
	@Test(dataProvider="src1")
	public void test01(final String inputFile,final String expr) 
		throws IOException
		{
		try {
		final Path output = support.createTmpPath(".vcf");

		
        Assert.assertEquals(new VCFStripAnnotations().instanceMain(new String[] {
        	"-o",output.toString(),
        	"-x",expr,
        	inputFile
        }),0);
        support.assertIsVcf(output);
		} finally {
			support.removeTmpFiles();
		}
		}
	
}

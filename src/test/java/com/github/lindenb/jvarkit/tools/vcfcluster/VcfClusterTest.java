package com.github.lindenb.jvarkit.tools.vcfcluster;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;

import htsjdk.samtools.util.IOUtil;

public class VcfClusterTest {

	@DataProvider(name = "src1")
	public Object[][] createData1() {
		final TestSupport support =new TestSupport();
		return new Object[][] {
			{support.resource("manta.B00I9CJ.vcf.gz"),10000},
			{support.resource("manta.B00GWGD.vcf.gz"),10000},
			{support.resource("manta.B00I9CJ.vcf.gz"),1000},
			{support.resource("manta.B00GWGD.vcf.gz"),1000},
			}
			;
		}
	
	@Test(dataProvider = "src1")
	public void testSV(final String invcf,int size) throws IOException {
		final TestSupport support =new TestSupport();
		Path dir = null;
		try {
			dir =  support.createTmpDirectory();
			
			Assert.assertEquals(new VcfCluster().instanceMain(new String[] {
				"-o",dir.toString(),
				"-n",String.valueOf(size),
				"--force",
				invcf
				}),0);
		} finally
			{
			if(dir!=null) IOUtil.deleteDirectoryTree(dir.toFile());
			support.removeTmpFiles();
			}
		}
}

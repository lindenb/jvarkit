package com.github.lindenb.jvarkit.tools.groupbygene;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;

public class GroupByGeneTest {
	
	private final TestSupport support = new TestSupport();
	
@DataProvider(name = "src1")
public Object[][] createData1() {
	return new Object[][]{
		{support.resource("ExAC.r1.sites.vep.vcf.gz")},
		{support.resource("gnomad.genomes.r2.0.1.sites.1.vcf.gz")},
		{support.resource("gnomad.exomes.r2.0.1.sites.vcf.gz")},
		{support.resource("rotavirus_rf.ann.vcf.gz")},
		};
	}

@Test(dataProvider="src1")
public void test01(final String vcf) throws IOException {
	try {
		final Path out = support.createTmpPath(".txt");
		
		Assert.assertEquals( new GroupByGene().instanceMain(new String[]{
				"-o",out.toString(),
				vcf
				}),0);
		Assert.assertTrue(support.wc(out)>1L);
		support.assertTsvTableIsConsitent(out, null);
		} 
	finally
		{
		support.removeTmpFiles();
		}
	}
}

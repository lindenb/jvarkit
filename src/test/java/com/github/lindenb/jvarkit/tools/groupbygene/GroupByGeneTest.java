package com.github.lindenb.jvarkit.tools.groupbygene;

import java.io.File;
import java.io.IOException;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

public class GroupByGeneTest extends TestUtils{
@DataProvider(name = "src1")
public Object[][] createData1() {
	return new Object[][]{
		{SRC_TEST_RESOURCE+"/ExAC.r1.sites.vep.vcf.gz"},
		{SRC_TEST_RESOURCE+"/gnomad.genomes.r2.0.1.sites.1.vcf.gz"},
		{SRC_TEST_RESOURCE+"/gnomad.exomes.r2.0.1.sites.vcf.gz"},
		{SRC_TEST_RESOURCE+"/rotavirus_rf.ann.vcf.gz"},
		};
	}

@Test(dataProvider="src1")
public void test01(final String vcf) throws IOException {
	final File out = super.createTmpFile(".txt");
	Assert.assertEquals(
			new GroupByGene().instanceMain(newCmd().
			add("-o").add(out).
			add(vcf).
			make()
			),0);
	Assert.assertTrue(wc(out)>1L);
	assertTsvTableIsConsitent(out, null);
	}
}

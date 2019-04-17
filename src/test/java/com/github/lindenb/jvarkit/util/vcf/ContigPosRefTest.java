package com.github.lindenb.jvarkit.util.vcf;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;


import htsjdk.variant.variantcontext.Allele;

public class ContigPosRefTest {

	@DataProvider(name = "src1")
	public Object[][] createData1() {
		return new Object[][]{
				{"chr1",123,"A"}
		};
	}
	@Test(dataProvider="src1")
	public void test1(String chrom,int pos,String allele) {
		ContigPosRef cp1 = new ContigPosRef(chrom,pos,Allele.create(allele,true));
		Assert.assertEquals(cp1.getContig(), chrom);
		Assert.assertEquals(cp1.getStart(), pos);
		Assert.assertEquals(cp1.getEnd(), pos);
		Assert.assertTrue(cp1.getReference().isReference());
		Assert.assertTrue(cp1.getReference().getDisplayString().equals(allele));
		ContigPosRef cp2 = new ContigPosRef(chrom,pos+1,Allele.create(allele,true));
		Assert.assertEquals(cp2.getStart(), pos+1);
		Assert.assertNotEquals(cp1, cp2);
		Assert.assertTrue(cp1.compareTo(cp2)<0);
	}
}

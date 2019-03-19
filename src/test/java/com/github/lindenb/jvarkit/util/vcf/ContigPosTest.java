package com.github.lindenb.jvarkit.util.vcf;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;


public class ContigPosTest{

	@DataProvider(name = "src1")
	public Object[][] createData1() {
		return new Object[][]{
				{"chr1:123","chr1",123}
		};
	}
	@Test(dataProvider="src1")
	public void test1(final String s1,String chrom,int pos) {
		ContigPos cp1 = new ContigPos(s1);
		Assert.assertEquals(cp1.getContig(), chrom);
		Assert.assertEquals(cp1.getStart(), pos);
		Assert.assertEquals(cp1.getEnd(), pos);
		ContigPos cp2 = new ContigPos(chrom,pos);
		Assert.assertEquals(cp1, cp2);
		Assert.assertEquals(0, cp1.compareTo(cp2));
	}
}

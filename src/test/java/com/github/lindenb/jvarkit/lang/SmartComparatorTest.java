package com.github.lindenb.jvarkit.lang;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

public class SmartComparatorTest  {
	@DataProvider(name="t0")
	public Object[][] getTestData0() {
		return new Object[][] {
			{"ABC","bcd"},
			{"T001","t02"},
			{"1t","00002T"},
			{"1","00009999999999999999999999999999999999999999999999999999999999"}
		};
	}
	@Test(dataProvider="t0")
	public void testLT(final String s1,final String s2)
		{
		Assert.assertTrue(new SmartComparator().compare(s1, s2)<0, s1+" < "+s2);
		}
	@Test(dataProvider="t0")
	public void testGE(final String s1,String s2)
		{
		Assert.assertTrue(new SmartComparator().compare(s2, s1)>=0, s1+" >0 "+s2);
		}
	
	@DataProvider(name="t1")
	public Object[][] getTestData1() {
		return new Object[][] {
			{"T00000000000001","T1"},
			{"00000000000001t","1T"}
		};
	}
	
	@Test(dataProvider="t1")
	public void testEQ(final String s1,final String s2)
		{
		Assert.assertTrue(new SmartComparator().compare(s1, s2) == 0, s1+" == "+s2);
		}
	
}

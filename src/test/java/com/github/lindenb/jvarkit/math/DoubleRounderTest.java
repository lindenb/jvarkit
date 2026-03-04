package com.github.lindenb.jvarkit.math;

import org.testng.Assert;
import org.testng.annotations.Test;

public class DoubleRounderTest {
	@Test
	public void test1() {
		DoubleRounder r = new DoubleRounder(2);
		double v2 =  r.applyAsDouble(0.1);
		Assert.assertTrue(v2>=0.1);
		Assert.assertTrue(v2<=0.101);
		}
}
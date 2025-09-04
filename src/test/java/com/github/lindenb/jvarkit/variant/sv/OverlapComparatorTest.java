package com.github.lindenb.jvarkit.variant.sv;

import org.testng.Assert;
import org.testng.annotations.Test;


public class  OverlapComparatorTest {
	@Test
	public void testOverlaps() {
		Assert.assertNotNull(OverlapComparator.makeDefault());
		
		OverlapComparator oc=OverlapComparator.parse("5:1%;100:50%;1000:99%");
		Assert.assertNotNull(oc);
		
		Assert.assertTrue(oc.test(100, 200, 101, 199));
		Assert.assertFalse(oc.test(100, 200, 201, 300));
		Assert.assertFalse(oc.test(100, 200, 199, 300));
		Assert.assertTrue(oc.test(1,1,1,2));
		Assert.assertFalse(oc.test(1, 1000000, 1, 100000));
		}
	
	
}

package com.github.lindenb.jvarkit.samtools.util;

import org.testng.Assert;
import org.testng.annotations.Test;

public class SimpleIntervalTest {
@Test
public void test01() {
	SimpleInterval r= new SimpleInterval("chr1:1-100");
	Assert.assertEquals(r.getContig(), "chr1");
	Assert.assertEquals(r.getStart(), 1);
	Assert.assertEquals(r.getEnd(), 100);
	Assert.assertEquals(r.length(), r.getLengthOnReference());
	
	SimpleInterval r2= new SimpleInterval(r);

	Assert.assertEquals(r,r2);
	}
}
